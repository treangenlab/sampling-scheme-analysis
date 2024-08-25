#!/usr/bin/env python

import os
import argparse
from fractions import Fraction
from models import *

# Needed for running on slurm
import sys
sys.path.append(os.getcwd()) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--window-size", nargs="+", required=True) 
    parser.add_argument("-k", "--kmer", nargs="+", required=True) 
    parser.add_argument("--sigma", nargs="+", required=True) 
    parser.add_argument("--local", help="Find minimum local scheme density", action="store_true") 
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-v", "--verbose", help="Log ILP to output", action="store_true")
    args = parser.parse_args()

    ws = args.window_size
    ks = args.kmer
    sigmas = args.sigma
    is_local = args.local
    n_threads = args.threads
    out_dir = "local" if is_local else "fwd"
    param_list = sorted(itertools.product(ws, ks, sigmas), key=lambda x: int(x[0]) + int(x[1]))

    for w, k, sigma in param_list:
        w = int(w)
        k = int(k)
        sigma = int(sigma)
        print(f"Working on w={w}, k={k}, s={sigma}")
        num_contexts = sigma ** (2*w + k - 1 if is_local else w + k)
        if num_contexts > 2**16:
            print("Too many contexts to feasibly solve...")
            continue

        # Set model
        gp.setParam("LogToConsole", args.verbose)
        wksig_to_sols = {}
        wksig_to_dens = {}
            
        gp.setParam("LogFile", f"{out_dir}/logs/w{w}-k{k}-s{sigma}.log")

        if not is_local:
            if os.path.isfile(f"fwd/sols/w{w}-k{k}-s{sigma}.pck"):
                print("Already completed!")
                # exit(0)

            seed = None
            alphabet = list(str(c) for c in range(sigma))
            for i in range(1, k-1):
                if os.path.isfile(f"fwd/sols/w{w}-k{k-i}-s{sigma}.pck"):
                    alphabet = list(str(c) for c in range(sigma))
                    with open(f"fwd/sols/w{w}-k{k-i}-s{sigma}.pck", 'rb') as pck_in:
                        seed = pck.load(pck_in)[(w, k-i, sigma)][0]
                    
                    suffices = list("".join(x) for x in itertools.product(alphabet, repeat=i))
                    seed_ext = {}
                    for node, val in seed.items():
                        seed_ext.update({node + c : val for c in suffices})
                    seed = seed_ext
                    break
                    
            if seed is None and os.path.isfile(f"fwd/sols/w{w-1}-k{k}-s{sigma}.pck"):
                with open(f"fwd/sols/w{w-1}-k{k}-s{sigma}.pck", 'rb') as pck_in:
                    seed = pck.load(pck_in)[(w-1, k, sigma)][0]
                    seed_ext = {}
                    for node, val in seed.items():
                        seed_ext.update({node + c : val for c in alphabet})
                    seed = seed_ext

            m = get_ILP_fwd(w, k, sigma=sigma, heuristics=1 if (k % w) == 1 else 0.9, seed=seed)

        else:
            if not os.path.isfile(f"fwd/sols/w{w}-k{k}-s{sigma}.pck"):
                print(f"({w}, {k}, {sigma})-forward scheme not found. There's no hope...")
                continue

            with open(f"fwd/sols/w{w}-k{k}-s{sigma}.pck", 'rb') as pck_in:
                fwd_seed = pck.load(pck_in)[w, k, sigma][0]

            with open(f"fwd/dens/w{w}-k{k}-s{sigma}.pck", 'rb') as pck_in:
                fwd_dens = pck.load(pck_in)[w, k, sigma]

            m = get_ILP_local(w, k, sigma=sigma, seed=fwd_seed)

        # Optimize model
        gp.setParam("Threads", n_threads)
        m.optimize()

        if m.status != GRB.OPTIMAL:
            print(f"Optimization terminated at non-zero gap = {100*m.MIPGap:.2f}%")
            continue
        
        sampled = int(np.round(m.ObjVal))
        density = Fraction(sampled, num_contexts)
        wksig_to_dens[(w, k, sigma)] = density
        if is_local and density < fwd_dens:
            print("({w}, {k}, {sigma})-local has lower density than fwd!")
        print(f"w={w}, k={k}, sigma={sigma}  {sampled}/{num_contexts}  -->  {float(density)}")

        dens_pck = f"{out_dir}/dens/w{w}-k{k}-s{sigma}.pck"
        with open(dens_pck, 'wb') as dens_out:
            pck.dump(wksig_to_dens, dens_out)

        # Save fwd solutions
        if not is_local:
            nFound = m.SolCount
            sols = []
            names = [x.VarName for x in m.getVars()]

            for sol in (range(nFound)):
                gp.setParam("SolutionNumber", sol)
                Ys = defaultdict(int)
                positions = set()
                wm = {}
                for name, val in zip(names, m.Xn):
                    if "x_" in name:
                        wm[name.split("_")[1]] = int(np.round(val))
                sols.append(wm)

            wksig_to_sols[(w, k, sigma)] = sols
            sols_pck = f"fwd/sols/w{w}-k{k}-s{sigma}.pck"
            with open(sols_pck, 'wb') as sols_out:
                pck.dump(wksig_to_sols, sols_out)

