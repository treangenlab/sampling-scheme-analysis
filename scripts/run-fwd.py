#!/usr/bin/env python

import os, sys
sys.path.append(os.getcwd()) 
import argparse
from fractions import Fraction
import gurobipy as gp
from gurobipy import GRB

from utils import *


def get_ILP_fwd(n, r, sigma=2, nSolutions=1, heuristics=0.1, seed=None):
    pgap = 0
    k = n - r + 1
    # print(n, k, r, 2**(n+1) * (1 - aperiodic_bound_suff(w, k)))

    # gp.setParam("BestObjStop", math.floor(sigma**(w + k)) * (1 - aperiodic_bound_suff(w, k, sigma)))
    gp.setParam("PoolSearchMode", 2)
    gp.setParam("PoolSolutions", nSolutions)
    gp.setParam("Heuristics", heuristics)
    gp.setParam("PoolGap", pgap)

    alphabet = list(str(c) for c in range(sigma))
    nodes = list("".join(x) for x in itertools.product(alphabet, repeat=n))
    edges = [(x, x[1:] + b) for x in nodes for b in alphabet]

    try:
        # Create a new model
        m = gp.Model("mip1")

        # Create variables
        x = {node: m.addVar(lb=0, ub=r-1, vtype=GRB.INTEGER, name=f"x_{node}") for node in nodes}
        y = {(u, v): m.addVar(vtype=GRB.BINARY, name=f"y_{u+v[-1]}") for u,v in edges}

        for (u, v) in edges:
            m.addConstr((y[(u,v)] == 0) >> (x[u] == x[v] + 1))
            m.addConstr((y[(u,v)] == 1) >> (x[u] <= x[v]))

        if seed:
            for node in nodes:
                x[node].Start = seed[node]


        
        if k % r == 1:
            for x in range(n+1, n+2):
                necklaces = get_necklaces(x, sigma)
                for neck, rots in necklaces.items():
                    p = len(rots)
                    neck_edges = [(rots[i][:n], rots[(i+1) % p][:n]) for i in range(p)]
                    m.addConstr(sum(y[e] for e in neck_edges) >= math.ceil(p / r))
        else:
            for d in range(1, min(2*n, 18), 1):
                necks = [neck for neck, rots in get_necklaces(d, sigma).items() if len(rots) & r == 1]
                for neck in necks:
                    s = ""
                    while len(s) < 10*(n+1):
                        s += neck
                    cycle = [s[:n+1]]
                    i = 1
                    while s[i:i+n+1] != cycle[0]:
                        cycle.append(s[i:i+n+1])
                        i += 1
                    
                    p = len(cycle)
                    neck_edges = [(v[:-1], v[1:]) for v in cycle]
                    m.addConstr(sum(y[e] for e in neck_edges) >= math.ceil(p / r))
                
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
        return

    except AttributeError as e:
        print('Encountered an attribute error')
        raise e

    # Set objective
    m.setObjective(sum(y.values()), GRB.MINIMIZE)
            
    return m


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--window-size", type=int, required=True) 
    parser.add_argument("-k", "--kmer", type=int, required=True) 
    parser.add_argument("--sigma", type=int, required=True) 
    parser.add_argument("-v", "--verbose", help="Log ILP to output", action="store_true")
    args = parser.parse_args()

    w = args.window_size
    k = args.kmer
    sigma = args.sigma

    # Set model
    gp.setParam("LogToConsole", args.verbose)
    wksig_to_sols = {}
    wksig_to_dens = {}
        
    n = w + k - 1
    gp.setParam("LogFile", f"fwd/logs/w{w}-k{k}-s{sigma}.log")

   
    if os.path.isfile(f"fwd/sols/w{w}-k{k}-s{sigma}.pck"):
        print("Already completed!")
        exit(0)

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
            suffices = list("".join(x) for x in itertools.product(alphabet, repeat=i))
            seed_ext = {}
            for node, val in seed.items():
                seed_ext.update({node + c : val for c in alphabet})
            seed = seed_ext

    m = get_ILP_fwd(n, w, sigma=sigma, heuristics=1 if (k % w) == 1 else 0.9, seed=seed)

    # Optimize model
    gp.setParam("Threads", 120)
    m.optimize()
    
    density = Fraction(int(np.round(m.ObjVal)), (sigma**(n+1)))
    wksig_to_dens[(w, k, sigma)] = density
    print(f"w={w}, k={k}, sigma={sigma}  {m.ObjVal}/{sigma**(n+1)}  -->  {float(density)}")
    
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

    dens_pck = f"fwd/dens/w{w}-k{k}-s{sigma}.pck"
    with open(dens_pck, 'wb') as dens_out:
        pck.dump(wksig_to_dens, dens_out)
