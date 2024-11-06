#!/usr/bin/env python

import os
import argparse
from fractions import Fraction
import re
from models import *

# Needed for running on slurm
import sys
sys.path.append(os.getcwd()) 


def get_solutions(model, is_local=False):
    nFound = model.SolCount
    sols = []
    names = [x.VarName for x in model.getVars()]

    for sol in (range(nFound)):
        gp.setParam("SolutionNumber", sol)
        Ys = defaultdict(int)
        positions = set()
        wm = {}
        for name, val in zip(names, model.Xn):
            if "x_" in name:
                if is_local:
                    if int(np.round(val)) == 1:
                        node, val = re.match(r'x_\((\d+), (\d+)\)', name).groups()
                        wm[node] = val
                else:
                    wm[name.split("_")[1]] = int(np.round(val))
        sols.append(wm)
    return sols

def get_fwd_solutions(model):
    nFound = model.SolCount
    sols = []
    names = [x.VarName for x in model.getVars()]

    for sol in (range(nFound)):
        gp.setParam("SolutionNumber", sol)
        Ys = defaultdict(int)
        positions = set()
        wm = {}
        for name, val in zip(names, model.Xn):
            if "x_" in name:
                wm[name.split("_")[1]] = int(np.round(val))
        sols.append(wm)
    return sols

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--window-size", nargs="+", required=True) 
    parser.add_argument("-k", "--kmer", nargs="+", required=True) 
    parser.add_argument("--sigma", nargs="+", required=True) 
    parser.add_argument("--local", help="Find minimum local scheme density", action="store_true") 
    parser.add_argument("-o", "--output", type=str, help="Path to output directory")
    parser.add_argument("--time-limit", type=int, help="Time limit (in seconds)")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-v", "--verbose", help="Log ILP to output", action="store_true")
    args = parser.parse_args()

    ws = list(map(int, args.window_size))
    ks = list(map(int, args.kmer))
    sigmas = list(map(int, args.sigma))
    is_local = args.local
    n_threads = args.threads
    if args.output:
        out_dir = args.output
    else:
        out_dir = "local" if is_local else "fwd"

    for sub in ("dens", "sols", "logs", "models"):
        os.makedirs(f"{out_dir}/{sub}", exist_ok=True)
    
    # Directories for results where there is gap in ILP bound
    os.makedirs(f"{out_dir}/dens/gap", exist_ok=True)
    os.makedirs(f"{out_dir}/sols/gap", exist_ok=True)

    param_list = sorted(itertools.product(ws, ks, sigmas), key=lambda x: int(x[0]) + int(x[1]))

    last_failed = False
    last_was_tight = False
    for w, k, sigma in param_list:
        if not is_local and last_failed:
            if last_was_tight:
                break

            if k % w != 1:
                continue

        last_was_tight = k % w == 1
        if os.path.isfile(f"{out_dir}/sols/w{w}-k{k}-s{sigma}.pck"):
            print(f"w{w}-k{k}-s{sigma} already completed!")
            # continue

        print(f"Working on w={w}, k={k}, s={sigma}")
        num_contexts = sigma ** (2*w + k - 2 if is_local else w + k)

        # Set model
        gp.setParam("LogToConsole", args.verbose)
        wksig_to_sols = {}
        wksig_to_dens = {}
            
        gp.setParam("LogFile", f"{out_dir}/logs/w{w}-k{k}-s{sigma}.log")

        if not is_local:

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

            model = get_ILP_fwd(w, k, sigma=sigma, heuristics=1 if ((k % w) == 1 or w == 2) else 0.8, seed=seed, method=3, concurrentMIP=max(1, n_threads//2), time_limit=args.time_limit)

        else:
            if not os.path.isfile(f"fwd/sols/w{w}-k{k}-s{sigma}.pck"):
                print(f"({w}, {k}, {sigma})-forward scheme not found. There's no hope...")
                continue

            with open(f"fwd/sols/w{w}-k{k}-s{sigma}.pck", 'rb') as pck_in:
                fwd_seed = pck.load(pck_in)[w, k, sigma][0]

            with open(f"fwd/dens/w{w}-k{k}-s{sigma}.pck", 'rb') as pck_in:
                fwd_dens = pck.load(pck_in)[w, k, sigma]

            model = get_ILP_local(w, k, sigma=sigma, seed=fwd_seed, method=3, concurrentMIP=max(1, n_threads//32), time_limit=args.time_limit)

        # Optimize model
        gp.setParam("Threads", n_threads)
        model.write(f"{out_dir}/models/w{w}-k{k}-s{sigma}.mps")
        model.optimize()

        if model.status != GRB.OPTIMAL:
            print(f"Optimization terminated at non-zero gap = {100*model.MIPGap:.2f}%")
            last_failed = True
        else:
            last_failed = False

        sampled = int(np.round(model.ObjVal))
        density = Fraction(sampled, num_contexts)
        wksig_to_dens[(w, k, sigma)] = density
        if is_local and density < fwd_dens:
            print(f"({w}, {k}, {sigma})-local has lower density than fwd!")
        if model.status == GRB.OPTIMAL:
            print(f"w={w}, k={k}, sigma={sigma}  {sampled}/{num_contexts}  -->  {float(density)}")

        if model.status == GRB.OPTIMAL:
            dens_pck = f"{out_dir}/dens/w{w}-k{k}-s{sigma}.pck"
            sols_pck = f"{out_dir}/sols/w{w}-k{k}-s{sigma}.pck"
        elif is_local and density < fwd_dens:
            dens_pck = f"{out_dir}/dens/gap/w{w}-k{k}-s{sigma}.pck"
            sols_pck = f"{out_dir}/sols/gap/w{w}-k{k}-s{sigma}.pck"

        with open(dens_pck, 'wb') as dens_out:
            pck.dump(wksig_to_dens, dens_out)

        # Save solutions
        sols = get_solutions(model, is_local)

        wksig_to_sols[(w, k, sigma)] = sols
        with open(sols_pck, 'wb') as sols_out:
            pck.dump(wksig_to_sols, sols_out)

