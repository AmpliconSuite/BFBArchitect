from pulp import LpMinimize, LpProblem, LpStatus, lpSum, LpVariable, PULP_CBC_CMD
import sys

def reconstruct_BFB_cbc(C, L, R, start, max_time=900, max_threads=8):
    model = LpProblem(name="BFB_reconstruction", sense=LpMinimize)
    
    T = round(max(sum(L) + sum(R) + 1, max(C)))
    # Initialize binary variables for consecutive-sequences
    segment_num = len(C)
    cs = {}
    for i in range(1, segment_num+1):
        for j in range(i, segment_num+1):
            for t in range(1, T+1):
                key1 = f'{i}_{j}_0_{t}'
                key2 = f'{i}_{j}_1_{t}'
                cs[key1] = LpVariable(name=f'cs_{key1}', cat="Binary")
                cs[key2] = LpVariable(name=f'cs_{key2}', cat="Binary")
    
    # Define the objective function
    objective = 0
    segment_cn_error = []
    left_fb_error, right_fb_error = [], []
    w1, w2, w3 = 1, 1, segment_num/2  # weights for different error types
    for k in range(1, segment_num+1):
        segment_cn_error.append(LpVariable(name=f'segment_error_{k}', lowBound=0))
        objective += w1 * segment_cn_error[k-1]
        left_fb_error.append(LpVariable(name=f'left_fb_error_{k}', lowBound=0))
        objective += w2 * left_fb_error[k-1]
        right_fb_error.append(LpVariable(name=f'right_fb_error_{k}', lowBound=0))
        objective += w2 * right_fb_error[k-1]
        # missing foldbacks
        right_fb = 0
        for i in range(1, k+1):
            for t in range(1, T):
                right_fb += cs[f'{i}_{k}_0_{t}']
        if R[k-1] == 0:
            missing_fb = LpVariable(name=f'missing_right_fb_{k}', cat="Binary")
            # model += (missing_fb >= right_fb, f'constraint_missing_right_{k}')
            model += (missing_fb <= right_fb, f'constraint_missing_right_{k}_lower_bound')
            model += (T * missing_fb >= right_fb, f'constraint_missing_right_{k}_upper_bound')
            objective += w3 * missing_fb
        left_fb = 0
        for j in range(k, segment_num+1):
            for t in range(1, T):
                left_fb += cs[f'{k}_{j}_1_{t}']
        if L[k-1] == 0:
            missing_fb = LpVariable(name=f'missing_left_fb_{k}', cat="Binary")
            # model += (missing_fb >= left_fb, f'constraint_missing_left_{k}')
            model += (missing_fb <= left_fb, f'constraint_missing_left_{k}_lower_bound')
            model += (T * missing_fb >= left_fb, f'constraint_missing_left_{k}_upper_bound')
            objective += w3 * missing_fb

    model += objective 

    for k in range(1, segment_num+1):
        # Segment copy number discrepancy
        seg_cn = 0
        for i in range(1, k+1):
            for j in range(k, segment_num+1):
                for t in range(1, T+1):
                    seg_cn += cs[f'{i}_{j}_0_{t}'] + cs[f'{i}_{j}_1_{t}']
        model += (seg_cn - C[k-1] <= segment_cn_error[k-1], f'constraint_segment_{k}')
        model += (seg_cn - C[k-1] >= -segment_cn_error[k-1], f'constraint_segment_{k}*')
        # Right foldback count discrepancy
        right_fb = 0
        for i in range(1, k+1):
            for t in range(1, T):
                right_fb += cs[f'{i}_{k}_0_{t}']
        model += (right_fb - R[k-1] <= right_fb_error[k-1], f'constraint_right_{k}')
        model += (right_fb - R[k-1] >= -right_fb_error[k-1], f'constraint_right_{k}*')
        # Left foldback count discrepancy
        left_fb = 0
        for j in range(k, segment_num+1):
            for t in range(1, T):
                left_fb += cs[f'{k}_{j}_1_{t}']
        model += (left_fb - L[k-1] <= left_fb_error[k-1], f'constraint_left_{k}')
        model += (left_fb - L[k-1] >= -left_fb_error[k-1], f'constraint_left_{k}*')
    
    # Add the constraints to the model
    # Define variables for constraints
    pa = {} # indicator for palindrome between time t1 and t2 (exclusive)
    rc = {} # indicator for reverse complement between consecutive-sequences added at t1 + 1 and t2 - 1
    sc = {} # indicator for super consecutive-sequence added at time t1 for consecutive-sequence added at time t2 and a palindrome between t1 and t2
    keys = [f'{i}_{j}_0' for i in range(1, segment_num+1) for j in range(i, segment_num+1)] + \
                    [f'{i}_{j}_1' for i in range(1, segment_num+1) for j in range(i, segment_num+1)]
    for t1 in range(1, T+1):
        for t2 in range(t1+1, T+1):
            pa[f'{t1}_{t2}'] = LpVariable(name=f'pa_{t1}_{t2}', cat="Binary")
            for s in keys:
                rc[f'{s}_{t1}_{t2}'] = LpVariable(name=f'rc_{s}_{t1}_{t2}', cat="Binary")
                sc[f'{s}_{t1}_{t2}'] = LpVariable(name=f'sc_{s}_{t1}_{t2}', cat="Binary")
    # Constraint 1: Initial sequence
    if start > 0:
        model += (cs[f'1_{segment_num}_0_1'] == 1, 'constraint_initial_+')
    else:
        model += (cs[f'1_{segment_num}_1_1'] == 1, 'constraint_initial_-')
    # Constraint 2: Empty sequence is palindromic
    for t in range(1, T):
        model += (pa[f'{t}_{t+1}'] == 1, f'constraint_empty_palindrome_{t}_{t+1}')
    # Constraint 3: Exactly one consecutive-sequence is added at each time point
    for t in range(1, T+1):
        expression = lpSum([cs[f'{s}_{t}'] for s in keys])
        model += (expression == 1, f'constraint_one_cs_{t}')
    # Constraint 4: The direction alternates between consecutive time points
    for t in range(1, T):
        expression = lpSum([cs[f'{i}_{j}_0_{t}'] - cs[f'{i}_{j}_1_{t}'] for i in range(1, segment_num+1) for j in range(i, segment_num+1)]) + \
                     lpSum([cs[f'{i}_{j}_0_{t+1}'] - cs[f'{i}_{j}_1_{t+1}'] for i in range(1, segment_num+1) for j in range(i, segment_num+1)])
        model += (expression == 0, f'constraint_alternate_direction_{t}')
    # Constraint 5: For time points t1 and t2 (t1 < t2), 
    # if the consecutive-sequence added at t1 + 1 and t2 - 1 are reverse complements, rc[t1][t2] = 1
    for t1 in range(1, T):
        for t2 in range(t1+1, T+1):
            for s in keys:
                s_bar = f'{s[:-1]}{"1" if s[-1]=="0" else "0"}' # reverse complement of s
                model += (rc[f'{s}_{t1}_{t2}'] <= cs[f'{s}_{t1 + 1}'], f'constraint_rc1_{s}_{t1}_{t2}')
                model += (rc[f'{s}_{t1}_{t2}'] <= cs[f'{s_bar}_{t2 - 1}'], f'constraint_rc2_{s}_{t1}_{t2}')
                model += (rc[f'{s}_{t1}_{t2}'] >= cs[f'{s}_{t1 + 1}'] + cs[f'{s_bar}_{t2 - 1}'] - 1, f'constraint_rc3_{s}_{t1}_{t2}')
    # Constraint 6: For time points t1 and t2 (t1 < t2 and t2 - t1 - 1 is even), 
    # if the sequence betwwen t1 and t2 is a plindrome 
    # and the consecutive-sequence added at t1 is a reverse complement of that added at t2
    # then pa[t1][t2] = 1
    for t1 in range(1, T):
        for t2 in range(t1+3, T+1, 2):
            M = lpSum([rc[f'{s}_{t1}_{t2}'] for s in keys])
            model += (pa[f'{t1}_{t2}'] <= pa[f'{t1 + 1}_{t2 - 1}'], f'constraint_pa1_{t1}_{t2}')
            model += (pa[f'{t1}_{t2}'] <= M, f'constraint_pa2_{t1}_{t2}')
            model += (pa[f'{t1}_{t2}'] >= pa[f'{t1 + 1}_{t2 - 1}'] + M - 1, f'constraint_pa3_{t1}_{t2}')
    # Constraint 7: For time points t1 and t2 (t1 < t2), 
    # if the consecutive-sequence added at t1 is a super consecutive-sequence of that added at t2
    # and the sequence between t1 and t2 is a palindrome, then sc[t1][t2] = 1
    for t1 in range(1, T):
        for t2 in range(t1+1, T+1):
            for s in keys:
                i, j, d = s.split('_')
                if d == '0':
                    super_keys = [f'{i}_{k}_1' for k in range(int(j), segment_num+1)]
                else:
                    super_keys = [f'{k}_{j}_0' for k in range(1, int(i)+1)]
                super_count = lpSum([cs[f'{sk}_{t1}'] for sk in super_keys])
                model += (sc[f'{s}_{t1}_{t2}'] <= super_count, f'constraint_sc1_{s}_{t1}_{t2}')
                model += (sc[f'{s}_{t1}_{t2}'] <= pa[f'{t1}_{t2}'], f'constraint_sc2_{s}_{t1}_{t2}')
                model += (sc[f'{s}_{t1}_{t2}'] >= super_count + pa[f'{t1}_{t2}'] - 1, f'constraint_sc3_{s}_{t1}_{t2}')
    # Constraint 8: cs is bounded by sc
    for t2 in range(2, T+1):
        for s in keys:
            total_sc = lpSum([sc[f'{s}_{t1}_{t2}'] for t1 in range(1, t2)])
            model += (cs[f'{s}_{t2}'] <= total_sc, f'constraint_cs_sc_{s}_{t2}')
    
    # Solve the problem
    status = model.solve(PULP_CBC_CMD(timeLimit=max_time, threads=max_threads, msg=0, options=["RandomS 42"]))
    # print(f'Objective: {model.objective.value()}')

    # Get consecutive-sequences added at each time point (all cs = 1)
    consecutive_sequences = []
    for var in model.variables():
        if var.value() > 0 and var.name.startswith('cs_'):
            # print(var)
            i, j, d, t = var.name[3:].split('_')
            consecutive_sequences.append((int(i), int(j), int(d), int(t)))
        # elif 'error' in var.name:
        #     print(var, var.value())
    consecutive_sequences = sorted(consecutive_sequences, key=lambda x: x[3])
    # print(consecutive_sequences)
    BFB_string = []
    for (i, j, d, _) in consecutive_sequences:
        if d == 0:
            segment = [x for x in range(i, j+1)]
        else:
            segment = [-x for x in range(j, i-1, -1)]
        BFB_string += segment
    return BFB_string, model.objective.value()

def reconstruct_BFB_gurobi(C, L, R, start, max_time=900, max_threads=8, pool_solutions=50, log_file=None, verbose=False):
    import gurobipy as gp
    from gurobipy import GRB
    # Build model with a configured environment so Gurobi messages go to log and/or stdout as requested
    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 1 if (log_file or verbose) else 0)
    env.setParam('LogToConsole', 1 if verbose else 0)
    if log_file:
        env.setParam('LogFile', log_file)
    env.start()
    m = gp.Model("BFB_reconstruction", env=env)
    m.Params.TimeLimit = max_time
    m.Params.Threads = max_threads

    # Solution pool settings
    # Mode 2 searches for multiple solutions; PoolSolutions caps how many are kept.
    m.Params.PoolSearchMode = 2
    m.Params.PoolSolutions = pool_solutions
    # If you want only optimal solutions, keep PoolGap at 0 (default)
    m.Params.PoolGap = 0.0

    segment_num = len(C)
    T = round(max(sum(L) + sum(R) + 1, max(C)))

    # Variables
    cs = {}  # consecutive-sequence binary vars: cs[i_j_d_t]
    for i in range(1, segment_num + 1):
        for j in range(i, segment_num + 1):
            for t in range(1, T + 1):
                key0 = f"{i}_{j}_0_{t}"
                key1 = f"{i}_{j}_1_{t}"
                cs[key0] = m.addVar(vtype=GRB.BINARY, name=f"cs_{key0}")
                cs[key1] = m.addVar(vtype=GRB.BINARY, name=f"cs_{key1}")

    # Objective vars and weights
    w1, w2, w3 = 1, 1, segment_num / 2
    segment_cn_error = [m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name=f"segment_error_{k}") for k in range(1, segment_num + 1)]
    left_fb_error = [m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name=f"left_fb_error_{k}") for k in range(1, segment_num + 1)]
    right_fb_error = [m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name=f"right_fb_error_{k}") for k in range(1, segment_num + 1)]

    # Missing foldback binaries when L[k-1] or R[k-1] == 0
    missing_right_fb = {}
    missing_left_fb = {}

    # Objective expression
    obj = gp.LinExpr()

    # Accumulate objective components and create missing fb vars/constraints
    for k in range(1, segment_num + 1):
        obj.addTerms(w1, segment_cn_error[k - 1])
        obj.addTerms(w2, left_fb_error[k - 1])
        obj.addTerms(w2, right_fb_error[k - 1])

        # right_fb (sum of cs[i_k_0_t] for i=1..k, t=1..T-1)
        right_fb_expr = gp.quicksum(cs[f"{i}_{k}_0_{t}"] for i in range(1, k + 1) for t in range(1, T))
        if R[k - 1] == 0:
            missing_right_fb[k] = m.addVar(vtype=GRB.BINARY, name=f"missing_right_fb_{k}")
            # missing_right_fb[k] ≤ right_fb
            m.addConstr(missing_right_fb[k] <= right_fb_expr, name=f"constraint_missing_right_{k}_lower_bound")
            # T * missing_right_fb[k] ≥ right_fb
            m.addConstr(T * missing_right_fb[k] >= right_fb_expr, name=f"constraint_missing_right_{k}_upper_bound")
            obj.addTerms(w3, missing_right_fb[k])

        # left_fb (sum of cs[k_j_1_t] for j=k..segment_num, t=1..T-1)
        left_fb_expr = gp.quicksum(cs[f"{k}_{j}_1_{t}"] for j in range(k, segment_num + 1) for t in range(1, T))
        if L[k - 1] == 0:
            missing_left_fb[k] = m.addVar(vtype=GRB.BINARY, name=f"missing_left_fb_{k}")
            m.addConstr(missing_left_fb[k] <= left_fb_expr, name=f"constraint_missing_left_{k}_lower_bound")
            m.addConstr(T * missing_left_fb[k] >= left_fb_expr, name=f"constraint_missing_left_{k}_upper_bound")
            obj.addTerms(w3, missing_left_fb[k])

    # Set objective (minimize)
    m.setObjective(obj, GRB.MINIMIZE)

    # Discrepancy constraints per segment k
    for k in range(1, segment_num + 1):
        # seg_cn = sum cs[i_j_0_t] + cs[i_j_1_t] for i=1..k, j=k..segment_num, t=1..T
        seg_cn_expr = gp.quicksum(cs[f"{i}_{j}_0_{t}"] + cs[f"{i}_{j}_1_{t}"]
                                  for i in range(1, k + 1)
                                  for j in range(k, segment_num + 1)
                                  for t in range(1, T + 1))
        m.addConstr(seg_cn_expr - C[k - 1] <= segment_cn_error[k - 1], name=f"constraint_segment_{k}")
        m.addConstr(seg_cn_expr - C[k - 1] >= -segment_cn_error[k - 1], name=f"constraint_segment_{k}*")

        # Right foldback count discrepancy
        right_fb_expr = gp.quicksum(cs[f"{i}_{k}_0_{t}"] for i in range(1, k + 1) for t in range(1, T))
        m.addConstr(right_fb_expr - R[k - 1] <= right_fb_error[k - 1], name=f"constraint_right_{k}")
        m.addConstr(right_fb_expr - R[k - 1] >= -right_fb_error[k - 1], name=f"constraint_right_{k}*")

        # Left foldback count discrepancy
        left_fb_expr = gp.quicksum(cs[f"{k}_{j}_1_{t}"] for j in range(k, segment_num + 1) for t in range(1, T))
        m.addConstr(left_fb_expr - L[k - 1] <= left_fb_error[k - 1], name=f"constraint_left_{k}")
        m.addConstr(left_fb_expr - L[k - 1] >= -left_fb_error[k - 1], name=f"constraint_left_{k}*")

    # Auxiliary binary variables: pa, rc, sc
    pa = {}  # pa[t1_t2]
    rc = {}  # rc[s_t1_t2]
    sc = {}  # sc[s_t1_t2]

    keys = [f"{i}_{j}_0" for i in range(1, segment_num + 1) for j in range(i, segment_num + 1)] + \
           [f"{i}_{j}_1" for i in range(1, segment_num + 1) for j in range(i, segment_num + 1)]

    for t1 in range(1, T + 1):
        for t2 in range(t1 + 1, T + 1):
            pa_key = f"{t1}_{t2}"
            pa[pa_key] = m.addVar(vtype=GRB.BINARY, name=f"pa_{pa_key}")
            for s in keys:
                rc_key = f"{s}_{t1}_{t2}"
                sc_key = f"{s}_{t1}_{t2}"
                rc[rc_key] = m.addVar(vtype=GRB.BINARY, name=f"rc_{rc_key}")
                sc[sc_key] = m.addVar(vtype=GRB.BINARY, name=f"sc_{sc_key}")

    # Constraint 1: Initial sequence
    if start > 0:
        m.addConstr(cs[f"1_{segment_num}_0_1"] == 1, name="constraint_initial_+")
    else:
        m.addConstr(cs[f"1_{segment_num}_1_1"] == 1, name="constraint_initial_-")

    # Constraint 2: Empty sequence is palindromic
    for t in range(1, T):
        m.addConstr(pa[f"{t}_{t + 1}"] == 1, name=f"constraint_empty_palindrome_{t}_{t + 1}")

    # Constraint 3: Exactly one consecutive-sequence is added at each time point
    for t in range(1, T + 1):
        expr = gp.quicksum(cs[f"{s}_{t}"] for s in keys)
        m.addConstr(expr == 1, name=f"constraint_one_cs_{t}")

    # Constraint 4: Direction alternates between consecutive time points
    for t in range(1, T):
        expr = gp.quicksum(cs[f"{i}_{j}_0_{t}"] - cs[f"{i}_{j}_1_{t}"]
                           for i in range(1, segment_num + 1)
                           for j in range(i, segment_num + 1)) + \
               gp.quicksum(cs[f"{i}_{j}_0_{t + 1}"] - cs[f"{i}_{j}_1_{t + 1}"]
                           for i in range(1, segment_num + 1)
                           for j in range(i, segment_num + 1))
        m.addConstr(expr == 0, name=f"constraint_alternate_direction_{t}")

    # Constraint 5: Reverse complements between t1+1 and t2-1
    for t1 in range(1, T):
        for t2 in range(t1 + 1, T + 1):
            for s in keys:
                s_bar = f"{s[:-1]}{'1' if s[-1] == '0' else '0'}"
                # rc[s_t1_t2] == cs[s_{t1+1}] ∧ cs[s_bar_{t2-1}] via three inequalities
                m.addConstr(rc[f"{s}_{t1}_{t2}"] <= cs[f"{s}_{t1 + 1}"], name=f"constraint_rc1_{s}_{t1}_{t2}")
                m.addConstr(rc[f"{s}_{t1}_{t2}"] <= cs[f"{s_bar}_{t2 - 1}"], name=f"constraint_rc2_{s}_{t1}_{t2}")
                m.addConstr(rc[f"{s}_{t1}_{t2}"] >= cs[f"{s}_{t1 + 1}"] + cs[f"{s_bar}_{t2 - 1}"] - 1, name=f"constraint_rc3_{s}_{t1}_{t2}")

    # Constraint 6: Palindrome propagation
    for t1 in range(1, T):
        for t2 in range(t1 + 3, T + 1, 2):
            M = gp.quicksum(rc[f"{s}_{t1}_{t2}"] for s in keys)
            m.addConstr(pa[f"{t1}_{t2}"] <= pa[f"{t1 + 1}_{t2 - 1}"], name=f"constraint_pa1_{t1}_{t2}")
            m.addConstr(pa[f"{t1}_{t2}"] <= M, name=f"constraint_pa2_{t1}_{t2}")
            m.addConstr(pa[f"{t1}_{t2}"] >= pa[f"{t1 + 1}_{t2 - 1}"] + M - 1, name=f"constraint_pa3_{t1}_{t2}")

    # Constraint 7: Super consecutive sequence and palindrome linkage
    for t1 in range(1, T):
        for t2 in range(t1 + 1, T + 1):
            for s in keys:
                i, j, d = s.split("_")
                i, j = int(i), int(j)
                if d == "0":
                    super_keys = [f"{i}_{k}_1" for k in range(j, segment_num + 1)]
                else:
                    super_keys = [f"{k}_{j}_0" for k in range(1, i + 1)]
                super_count = gp.quicksum(cs[f"{sk}_{t1}"] for sk in super_keys)
                m.addConstr(sc[f"{s}_{t1}_{t2}"] <= super_count, name=f"constraint_sc1_{s}_{t1}_{t2}")
                m.addConstr(sc[f"{s}_{t1}_{t2}"] <= pa[f"{t1}_{t2}"], name=f"constraint_sc2_{s}_{t1}_{t2}")
                m.addConstr(sc[f"{s}_{t1}_{t2}"] >= super_count + pa[f"{t1}_{t2}"] - 1, name=f"constraint_sc3_{s}_{t1}_{t2}")

    # Constraint 8: cs at t2 must be supported by some sc from earlier t1
    for t2 in range(2, T + 1):
        for s in keys:
            total_sc = gp.quicksum(sc[f"{s}_{t1}_{t2}"] for t1 in range(1, t2))
            m.addConstr(cs[f"{s}_{t2}"] <= total_sc, name=f"constraint_cs_sc_{s}_{t2}")

    # Optimize and collect solution pool
    m.optimize()

    if m.Status not in [GRB.OPTIMAL, GRB.INTERRUPTED, GRB.TIME_LIMIT]:
        # Return empty if infeasible/unbounded etc.
        return [], None

    # Best objective value in the pool
    obj_value = m.ObjVal

    # Iterate through solutions in the pool
    n_solutions = m.SolCount
    solutions = []

    # Helper to build BFB string from a solution number
    def extract_BFB_string(solution_number):
        # Switch to specific solution in the pool
        m.Params.SolutionNumber = solution_number
        # Collect cs vars equal to 1 in this solution
        consecutive_sequences = []
        for v in m.getVars():
            if v.VarName.startswith("cs_"):
                # Xn is the value in the nth solution of the pool
                if v.Xn > 0.5:
                    # Parse name cs_i_j_d_t
                    _, payload = v.VarName.split("cs_")
                    i, j, d, t = payload.split("_")
                    consecutive_sequences.append((int(i), int(j), int(d), int(t)))
        # Sort by time t
        consecutive_sequences.sort(key=lambda x: x[3])
        # Build BFB string
        BFB_string = []
        for (i, j, d, _) in consecutive_sequences:
            if d == 0:
                segment = list(range(i, j + 1))
            else:
                segment = [-x for x in range(j, i - 1, -1)]
            BFB_string += segment
        return BFB_string

    # Collect only solutions with optimal objective value
    # PoolObjVal gives the objective of the current solution in the pool
    BFB_strings = set()
    for si in range(n_solutions):
        m.Params.SolutionNumber = si
        pool_obj = m.PoolObjVal
        # Keep only optimal solutions (exact match)
        if abs(pool_obj - obj_value) < 1e-9:
            BFB_list = extract_BFB_string(si)
            BFB_string = print_BFB_string(BFB_list, print_to_console=False)
            if BFB_string not in BFB_strings:
                solutions.append(BFB_list)
                BFB_strings.add(BFB_string)

    return solutions, obj_value

def reconstruct_BFB_mosek(C, L, R, start, max_time=900, max_threads=8, log_file=None, verbose=False):
    import mosek
    from mosek.fusion import Model, Domain, Expr, ObjectiveSense
    """
    Reconstruct BFB using MOSEK optimization solver.
    
    Args:
        C: Copy number vector
        L: Left foldback counts
        R: Right foldback counts  
        start: Initial orientation (>0 for +, <=0 for -)
        max_time: Time limit in seconds
        max_threads: Number of threads
        log_file: Path to log file
        verbose: Whether to print solver output
        
    Returns:
        solutions: List containing single optimal BFB string
        obj_value: Optimal objective value
    """
    
    m = Model("BFB_reconstruction")
        
    # Configure logging
    if verbose:
        m.setLogHandler(sys.stdout)
    else:
        m.setLogHandler(None)
    if log_file:
        m.setLogHandler(open(log_file, 'a'))
        
    
    # Set time limit (in seconds)
    m.setSolverParam("optimizerMaxTime", float(max_time))
    
    # Set number of threads
    m.setSolverParam("numThreads", max_threads)
    
    segment_num = len(C)
    T = round(max(sum(L) + sum(R) + 1, max(C)))
    
    # ============================================================
    # Variables
    # ============================================================
    
    # Binary variables cs[i,j,d,t]
    cs = {}
    for i in range(1, segment_num + 1):
        for j in range(i, segment_num + 1):
            for d in [0, 1]:
                for t in range(1, T + 1):
                    key = f"{i}_{j}_{d}_{t}"
                    cs[key] = m.variable(f"cs_{key}", Domain.binary())
    
    # Continuous error variables
    segment_cn_error = m.variable("segment_cn_error", segment_num, Domain.greaterThan(0.0))
    left_fb_error = m.variable("left_fb_error", segment_num, Domain.greaterThan(0.0))
    right_fb_error = m.variable("right_fb_error", segment_num, Domain.greaterThan(0.0))
    
    # Missing foldback binaries
    missing_right_fb = {}
    missing_left_fb = {}
    
    # ============================================================
    # Objective Function
    # ============================================================
    
    w1, w2, w3 = 1.0, 1.0, segment_num / 2.0
    obj_terms = []
    
    # Main error terms
    obj_terms.append(Expr.mul(w1, Expr.sum(segment_cn_error)))
    obj_terms.append(Expr.mul(w2, Expr.sum(left_fb_error)))
    obj_terms.append(Expr.mul(w2, Expr.sum(right_fb_error)))
    
    # Missing foldback penalties
    for k in range(1, segment_num + 1):
        # Right foldback
        right_fb_vars = [cs[f"{i}_{k}_0_{t}"] for i in range(1, k + 1) for t in range(1, T)]
        right_fb_sum = Expr.sum(Expr.vstack(right_fb_vars)) if right_fb_vars else Expr.constTerm(0.0)
        
        if R[k - 1] == 0:
            missing_right_fb[k] = m.variable(f"missing_right_fb_{k}", Domain.binary())
            m.constraint(f"missing_right_{k}_lb", Expr.sub(missing_right_fb[k], right_fb_sum), Domain.lessThan(0.0))
            m.constraint(f"missing_right_{k}_ub", Expr.sub(Expr.mul(T, missing_right_fb[k]), right_fb_sum), Domain.greaterThan(0.0))
            obj_terms.append(Expr.mul(w3, missing_right_fb[k]))
        
        # Left foldback
        left_fb_vars = [cs[f"{k}_{j}_1_{t}"] for j in range(k, segment_num + 1) for t in range(1, T)]
        left_fb_sum = Expr.sum(Expr.vstack(left_fb_vars)) if left_fb_vars else Expr.constTerm(0.0)
        
        if L[k - 1] == 0:
            missing_left_fb[k] = m.variable(f"missing_left_fb_{k}", Domain.binary())
            m.constraint(f"missing_left_{k}_lb", Expr.sub(missing_left_fb[k], left_fb_sum), Domain.lessThan(0.0))
            m.constraint(f"missing_left_{k}_ub", Expr.sub(Expr.mul(T, missing_left_fb[k]), left_fb_sum), Domain.greaterThan(0.0))
            obj_terms.append(Expr.mul(w3, missing_left_fb[k]))
    
    # Set objective
    m.objective("obj", ObjectiveSense.Minimize, Expr.add(obj_terms))
    
    # ============================================================
    # Segment Discrepancy Constraints
    # ============================================================
    
    for k in range(1, segment_num + 1):
        # Segment copy number
        seg_cn_vars = [cs[f"{i}_{j}_{d}_{t}"] 
                        for i in range(1, k + 1) 
                        for j in range(k, segment_num + 1)
                        for d in [0, 1]
                        for t in range(1, T + 1)]
        seg_cn_sum = Expr.sum(Expr.vstack(seg_cn_vars))
        
        # |seg_cn_sum - C[k-1]| <= error
        m.constraint(f"seg_cn_{k}_pos", Expr.sub(Expr.sub(seg_cn_sum, C[k-1]), segment_cn_error.index(k-1)), Domain.lessThan(0.0))
        m.constraint(f"seg_cn_{k}_neg", Expr.add(Expr.sub(seg_cn_sum, C[k-1]), segment_cn_error.index(k-1)), Domain.greaterThan(0.0))
        
        # Right foldback
        right_fb_vars = [cs[f"{i}_{k}_0_{t}"] for i in range(1, k + 1) for t in range(1, T)]
        right_fb_sum = Expr.sum(Expr.vstack(right_fb_vars)) if right_fb_vars else Expr.constTerm(0.0)
        
        m.constraint(f"right_fb_{k}_pos", Expr.sub(Expr.sub(right_fb_sum, R[k-1]), right_fb_error.index(k-1)), Domain.lessThan(0.0))
        m.constraint(f"right_fb_{k}_neg", Expr.add(Expr.sub(right_fb_sum, R[k-1]), right_fb_error.index(k-1)), Domain.greaterThan(0.0))
        
        # Left foldback
        left_fb_vars = [cs[f"{k}_{j}_1_{t}"] for j in range(k, segment_num + 1) for t in range(1, T)]
        left_fb_sum = Expr.sum(Expr.vstack(left_fb_vars)) if left_fb_vars else Expr.constTerm(0.0)
        
        m.constraint(f"left_fb_{k}_pos", Expr.sub(Expr.sub(left_fb_sum, L[k-1]), left_fb_error.index(k-1)), Domain.lessThan(0.0))
        m.constraint(f"left_fb_{k}_neg", Expr.add(Expr.sub(left_fb_sum, L[k-1]), left_fb_error.index(k-1)), Domain.greaterThan(0.0))
    
    # ============================================================
    # Auxiliary Variables
    # ============================================================
    
    pa = {}
    rc = {}
    sc = {}
    
    keys = [f"{i}_{j}_{d}" for i in range(1, segment_num + 1) 
            for j in range(i, segment_num + 1) for d in [0, 1]]
    
    for t1 in range(1, T + 1):
        for t2 in range(t1 + 1, T + 1):
            pa[f"{t1}_{t2}"] = m.variable(f"pa_{t1}_{t2}", Domain.binary())
            for s in keys:
                rc[f"{s}_{t1}_{t2}"] = m.variable(f"rc_{s}_{t1}_{t2}", Domain.binary())
                sc[f"{s}_{t1}_{t2}"] = m.variable(f"sc_{s}_{t1}_{t2}", Domain.binary())
    
    # ============================================================
    # BFB Constraints
    # ============================================================
    
    # Constraint 1: Initial sequence
    if start > 0:
        m.constraint("initial", cs[f"1_{segment_num}_0_1"], Domain.equalsTo(1.0))
    else:
        m.constraint("initial", cs[f"1_{segment_num}_1_1"], Domain.equalsTo(1.0))
    
    # Constraint 2: Empty palindrome
    for t in range(1, T):
        m.constraint(f"empty_pal_{t}", pa[f"{t}_{t+1}"], Domain.equalsTo(1.0))
    
    # Constraint 3: Exactly one cs per time
    for t in range(1, T + 1):
        cs_t = [cs[f"{s}_{t}"] for s in keys]
        m.constraint(f"one_cs_{t}", Expr.sum(Expr.vstack(cs_t)), Domain.equalsTo(1.0))
    
    # Constraint 4: Direction alternates
    for t in range(1, T):
        sum_t = Expr.sum(Expr.vstack([Expr.sub(cs[f"{i}_{j}_0_{t}"], cs[f"{i}_{j}_1_{t}"]) 
                            for i in range(1, segment_num + 1) for j in range(i, segment_num + 1)]))
        sum_t1 = Expr.sum(Expr.vstack([Expr.sub(cs[f"{i}_{j}_0_{t+1}"], cs[f"{i}_{j}_1_{t+1}"]) 
                            for i in range(1, segment_num + 1) for j in range(i, segment_num + 1)]))
        m.constraint(f"alternate_{t}", Expr.add(sum_t, sum_t1), Domain.equalsTo(0.0))
    
    # Constraint 5: Reverse complements
    for t1 in range(1, T):
        for t2 in range(t1 + 1, T + 1):
            for s in keys:
                parts = s.split("_")
                s_bar = f"{parts[0]}_{parts[1]}_{'1' if parts[2] == '0' else '0'}"
                
                rc_var = rc[f"{s}_{t1}_{t2}"]
                cs_t1 = cs[f"{s}_{t1+1}"]
                cs_t2 = cs[f"{s_bar}_{t2-1}"]
                
                # rc = cs_t1 AND cs_t2
                m.constraint(f"rc1_{s}_{t1}_{t2}", Expr.sub(rc_var, cs_t1), Domain.lessThan(0.0))
                m.constraint(f"rc2_{s}_{t1}_{t2}", Expr.sub(rc_var, cs_t2), Domain.lessThan(0.0))
                m.constraint(f"rc3_{s}_{t1}_{t2}", Expr.sub(Expr.add(cs_t1, cs_t2), Expr.add(rc_var, 1.0)), Domain.lessThan(0.0))
    
    # Constraint 6: Palindrome propagation
    for t1 in range(1, T):
        for t2 in range(t1 + 3, T + 1, 2):
            M = Expr.sum(Expr.vstack([rc[f"{s}_{t1}_{t2}"] for s in keys]))
            pa_curr = pa[f"{t1}_{t2}"]
            pa_inner = pa[f"{t1+1}_{t2-1}"]
            
            # pa_curr = pa_inner AND M
            m.constraint(f"pa1_{t1}_{t2}", Expr.sub(pa_curr, pa_inner), Domain.lessThan(0.0))
            m.constraint(f"pa2_{t1}_{t2}", Expr.sub(pa_curr, M), Domain.lessThan(0.0))
            m.constraint(f"pa3_{t1}_{t2}", Expr.sub(Expr.add(pa_inner, M), Expr.add(pa_curr, 1.0)), Domain.lessThan(0.0))
    
    # Constraint 7: Super consecutive
    for t1 in range(1, T):
        for t2 in range(t1 + 1, T + 1):
            for s in keys:
                i, j, d = map(int, s.split("_"))
                
                if d == 0:
                    super_keys = [f"{i}_{k}_1_{t1}" for k in range(j, segment_num + 1)]
                else:
                    super_keys = [f"{k}_{j}_0_{t1}" for k in range(1, i + 1)]
                
                super_sum = Expr.sum(Expr.vstack([cs[sk] for sk in super_keys]) if super_keys else Expr.constTerm(0.0))
                sc_var = sc[f"{s}_{t1}_{t2}"]
                pa_var = pa[f"{t1}_{t2}"]
                
                # sc = super_sum AND pa
                m.constraint(f"sc1_{s}_{t1}_{t2}", Expr.sub(sc_var, super_sum), Domain.lessThan(0.0))
                m.constraint(f"sc2_{s}_{t1}_{t2}", Expr.sub(sc_var, pa_var), Domain.lessThan(0.0))
                m.constraint(f"sc3_{s}_{t1}_{t2}", Expr.sub(Expr.add(super_sum, pa_var), Expr.add(sc_var, 1.0)), Domain.lessThan(0.0))
    
    # Constraint 8: cs must be supported by sc
    for t2 in range(2, T + 1):
        for s in keys:
            sc_sum = Expr.sum(Expr.vstack([sc[f"{s}_{t1}_{t2}"] for t1 in range(1, t2)]))
            m.constraint(f"cs_sc_{s}_{t2}", Expr.sub(cs[f"{s}_{t2}"], sc_sum), Domain.lessThan(0.0))
    
    # ============================================================
    # Solve
    # ============================================================
    
    m.solve()
    obj_value = m.primalObjValue()
    
    # ============================================================
    # Extract Solution
    # ============================================================
    
    consecutive_sequences = []
    for key, var in cs.items():
        if var.level()[0] > 0.5:
            parts = key.split("_")
            i, j, d, t = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
            consecutive_sequences.append((i, j, d, t))
    
    consecutive_sequences.sort(key=lambda x: x[3])
    
    BFB_string = []
    for (i, j, d, _) in consecutive_sequences:
        if d == 0:
            segment = list(range(i, j + 1))
        else:
            segment = [-x for x in range(j, i - 1, -1)]
        BFB_string += segment
    
    return [BFB_string], obj_value

def check_BFB_string(string):
    if len(string) == 0:
        return False
    original_len = max([abs(seg) for seg in string])
    pos = len(string) - 1
    # print(pos)
    while pos >= original_len:
        # Find the longest palindromic suffix
        positions = [i for i, seg in enumerate(string[:pos]) if seg == -string[pos]]
        i = None
        for j in positions:
            if is_palindrome(string[j:pos+1]):
                i = j
                break
        if i == None:
            return False
        pos = (i+pos)//2
        # print(pos)
    return True

def is_palindrome(sequence):
    # Check if a new sequence is a palindrome with even length
    if len(sequence) % 2 == 1:
        return False
    i, mid = 0, len(sequence)//2
    while i < mid:
        if sequence[i] != -sequence[-i-1]:
            return False
        i += 1
    return True

def print_BFB_string(BFB_string, print_to_console=True):
    sequence = []
    for segment in BFB_string:
        if segment > 0:
            sequence.append(str(segment) + '+')
        else:
            # sequence += str(-segment) + '\u0305 '
            sequence.append(str(-segment) + '-')
    string = ','.join(sequence)
    if print_to_console:
        print(string)
    return string

if __name__ == '__main__':
    # error-free examples
    start = 1
    C = [9, 5, 3, 4]
    L = [4, 0, 0, 1]
    R = [2, 1, 0, 2]

    start = -5
    C = [4, 4, 6, 3, 1]
    L = [2, 0, 2, 0, 0]
    R = [0, 1, 1, 1, 0]

    start = 1
    C = [1, 1, 5, 6, 8]
    L = [0, 0, 2, 0, 1]
    R = [0, 0, 0, 0, 4]

    start = 1
    C = [1, 7, 5]
    L = [0, 3, 0]
    R = [0, 1, 2]
    
    start = 1
    C = [1, 6, 8, 2]
    L = [0, 2, 1, 0]
    R = [0, 0, 3, 1]

    # examples with errors
    # start = 1
    # C = [1, 7, 5, 3]
    # L = [0, 3, 0, 0]
    # R = [0, 1, 1, 1]

    start = 1
    # C = [8, 6, 14, 34, 26, 23]
    # L = [3, 0, 0, 10, 0, 1]
    # R = [0, 0, 0, 5, 0, 12]

    # C = [2, 6, 11, 13, 9, 3]
    # L = [0, 2, 3, 0, 0, 0]
    # R = [0, 0, 0, 1, 3, 1]

    # C = [11, 9, 11, 5]
    # L = [4, 0, 0, 0]
    # R = [0, 0, 3, 1]

    # C = [11, 20, 22, 5]
    # L = [5, 5, 0, 0]
    # R = [0, 0, 6, 1]

    # C = [37, 34, 37, 41, 38, 39, 37, 38, 42, 40, 19, 10, 8, 11, 7]
    # L = [18, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    # R = [0, 0, 0, 0, 0, 1, 0, 0, 0, 12, 4, 0, 0, 2, 0]

    C = [0, 2, 2, 0, 2, 1]
    L = [0, 0, 0, 0, 0, 0]
    R = [0, 0, 1, 0, 0, 0]

    start = -5
    C = [5, 8, 9, 10, 9] 
    L = [0, 0, 0, 0, 0]
    R = [0, 0, 0, 0, 1]

    start = -8
    C = [6, 7, 6, 7, 11, 13, 9, 3]
    L = [2, 0, 0, 0, 3, 0, 0, 0]
    R = [0, 0, 0, 0, 0, 1, 3, 1]

    # BFB_string, obj_val = reconstruct_BFB_cbc(C, L, R, start, max_time=None, max_threads=12)
    # print('ILP objective value:', obj_val)
    # print('ILP output string', BFB_string)
    # print('BFB' if check_BFB_string(BFB_string) else 'Not BFB')
    # print_BFB_string(BFB_string)

    BFB_strings, obj_val = reconstruct_BFB_gurobi(C, L, R, start, max_time=900, max_threads=12, pool_solutions=50)
    print('Gurobi objective value:', obj_val)
    for idx, BFB_string in enumerate(BFB_strings):
        print(f'Gurobi output string {idx+1}:', BFB_string)
        print('BFB' if check_BFB_string(BFB_string) else 'Not BFB')
        print_BFB_string(BFB_string)