import gurobipy as gp
from gurobipy import GRB

from utils import *
TIME_LIMIT = 30 * 60 # Run for 30 min


def get_ILP_local(w, k, sigma=2, nSolutions=1, start=None, sketch_size=1):
    pgap = 0.0000
    context_size = 2*w + k - 1
    ell = w + k - 1
    
    gp.setParam("PoolSearchMode", 2)
    gp.setParam("PoolSolutions", nSolutions)
    gp.setParam("PoolGap", pgap)
    gp.setParam("Method", 3)
    gp.setParam("TimeLimit", TIME_LIMIT)

    
    alphabet = [str(c) for c in range(sigma)]
    nodes = list("".join(x) for x in itertools.product(alphabet, repeat=ell))
    edges = [(x, x[1:] + b) for x in nodes for b in alphabet]

    context_dbseq = de_bruijn_invBWT(context_size, alphabet)
    N = sigma ** context_size
    assert(N == len(context_dbseq))
    # allow for wraparound
    context_dbseq += context_dbseq
    
    try:
        # Create a new model
        m = gp.Model("mip1")
    
        # Create variables
        x = {node: [m.addVar(vtype=GRB.BINARY, name=f"x_({node}, {j})") for j in range(w)]
            for node in nodes}
        
        y = [m.addVar(vtype=GRB.BINARY , name=f"y_{a}") for a in range(N)]

        # Seed start values if provided
        if start:
            for node in x:
                for i in range(w):
                    x[node][i].Start = 0 if i != start[node] else 1
                
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    
    except AttributeError as e:
        print('Encountered an attribute error')
        raise e
        
    # Set objective
    m.setObjective(np.sum(y), GRB.MINIMIZE)

    # Ensure exactly sketch_size positions selected per window
    for x_i in x.values():
        m.addConstr(sum(x_i) == sketch_size)
        
    for pos in (range(N)):
        y_i = y[pos]
        # y_i is selected if any context selects it
        for j in range(0, w):
            W_pos = pos - j
            if W_pos < 0:
                W_pos += N
            W = context_dbseq[W_pos:W_pos+ell]
            m.addConstr(y_i >= x[W][j]) 
    return m


def get_ILP_fwd(w, k, sigma=2, nSolutions=1, heuristics=0.1, seed=None):
    pgap = 0
    gp.setParam("PoolSearchMode", 2)
    gp.setParam("PoolSolutions", nSolutions)
    gp.setParam("Heuristics", heuristics)
    gp.setParam("PoolGap", pgap)
    gp.setParam("Method", 3)
    gp.setParam("TimeLimit", TIME_LIMIT)

    n = w + k - 1
    alphabet = list(str(c) for c in range(sigma))
    nodes = list("".join(x) for x in itertools.product(alphabet, repeat=n))
    edges = [(x, x[1:] + b) for x in nodes for b in alphabet]

    try:
        # Create a new model
        m = gp.Model("mip1")

        # Create variables
        x = {node: m.addVar(lb=0, ub=w-1, vtype=GRB.INTEGER, name=f"x_{node}") for node in nodes}
        y = {(u, v): m.addVar(vtype=GRB.BINARY, name=f"y_{u+v[-1]}") for u,v in edges}

        for (u, v) in edges:
            m.addConstr((y[(u,v)] == 0) >> (x[u] == x[v] + 1))
            m.addConstr((y[(u,v)] == 1) >> (x[u] <= x[v]))

        if seed:
            for node in nodes:
                x[node].Start = seed[node]
        
        if k % w == 1:
            for x in range(n+1, n+2):
                necklaces = get_necklaces(x, sigma)
                for neck, rots in necklaces.items():
                    p = len(rots)
                    neck_edges = [(rots[i][:n], rots[(i+1) % p][:n]) for i in range(p)]
                    m.addConstr(sum(y[e] for e in neck_edges) >= math.ceil(p / w))
        else:
            cycle_constraints = 0

            # Reduce cycle space for larger alphabets
            max_cycle = {2: 16, 3: 12, 4: 8}[sigma]

            # Use other simple cycles
            for d in range(1, min(2*n, max_cycle) + 1):
                necks = [neck for neck, rots in get_necklaces(d, sigma).items()]
                for neck in necks:
                    s = ""
                    while len(s) < n+1:
                        s += neck
                    s = s + s
                    cycle = [s[:n+1]]
                    i = 1
                    while s[i:i+n+1] != cycle[0]:
                        cycle.append(s[i:i+n+1])
                        i += 1
                    
                    p = len(cycle)

                    # If not simple cycle
                    if p != len(set(cycle)):
                        continue

                    # If not helpful cycle
                    if p % w == 0:
                        continue

                    neck_edges = [(v[:-1], v[1:]) for v in cycle]
                    m.addConstr(sum(y[e] for e in neck_edges) >= math.ceil(p / w))
                    cycle_constraints += 1
            print(f"{cycle_constraints} cycle constraints added...")
                
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
        return

    except AttributeError as e:
        print('Encountered an attribute error')
        raise e

    # Set objective
    m.setObjective(sum(y.values()), GRB.MINIMIZE)
            
    return m

