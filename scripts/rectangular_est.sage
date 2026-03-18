# SPDX-License-Identifier: MIT
#
# Script to estimate SNOVA security against known attacks
#
# Copyright (c) 2026 SNOVA TEAM


import math
import sys
import traceback


variants = [
    [28, 5, 19, 4],
    [24, 5, 16, 4],
    [48, 17, 16, 2],
    [48, 16, 19, 2],
    [28, 4, 16, 4, 5],
    [28, 4, 19, 4, 5],

    [40, 7, 19, 4],
    [37, 8, 16, 4],
    [72, 25, 16, 2],
    [72, 24, 19, 2],
    [38, 5, 16, 4, 5],
    [38, 5, 19, 4, 5],

    [50, 9, 19, 4],
    [60, 10, 16, 4],
    [97, 33, 16, 2],
    [96, 32, 19, 2],
    [52, 6, 16, 4, 6],
    [52, 6, 19, 4, 6],
]


# SL I
variants1 = [
    [43, 6, 16, 2, 6, 16],
    [44, 5, 16, 2, 8, 16],
    [36, 5, 16, 3, 5],
    [26, 3, 16, 4, 6],
    [28, 4, 16, 4, 5, 5, 25],

    [43, 6, 19, 2, 6, 16],
    [43, 5, 23, 2, 7, 16],
    [43, 5, 19, 2, 8, 16],
    [28, 4, 19, 4, 5, 5, 25],
]

# SL3
variants3 = [
    [64, 8, 16, 2, 6, 24, 24],
    [48, 8, 16, 3, 4],
    [38, 5, 16, 4, 5],
    [44, 6, 16, 4, 5],

    [64, 8, 19, 2, 6],
    [38, 4, 19, 4, 6],
    [38, 5, 19, 4, 5],
    [44, 6, 19, 4, 5],
]

# SL 5
variants5 = [
    [86, 11, 16, 2, 6, 32, 30],
    [64, 9, 16, 3, 5],
    [52, 6, 16, 4, 6],
    [52, 7, 16, 4, 5],
    [56, 8, 16, 4, 5],

    [52, 6, 19, 4, 6],
    [52, 7, 19, 4, 5],
    [56, 8, 19, 4, 5],
]

variants += variants1 + variants3 + variants5


def RegDim(np, mp, q):
    precision = max(120, np + 1)
    R = PowerSeriesRing(ZZ, 't0', default_prec=precision)
    t0 = R.gens()[0]

    if True:
        pol = ((1 - t0**2) / (1 - t0**(2 * q)))**mp * ((1 - t0**q) / (1 - t0))**np / (1 - t0)
    else:
        pol = (1 - t0**2)**mp / (1 - t0)**(np + 1)

    pol_coef = pol.list()
    for d_reg in range(len(pol_coef)):
        if pol_coef[d_reg] <= 0:
            return d_reg
    return 9999


def MQ(n, m, q):
    dim = RegDim(n, m, q)
    if dim < 999:
        est = 3 * math.comb(n + 2, 2) * math.comb(n + dim, dim)**2
        return est
    else:
        return 2**4096


def HybridMQ(n, m, q):
    # Overdetermined case
    min = MQ(m, m, q)
    for k in range(max(0, n - m), n):
        try:
            est = MQ(n - k, m, q) * q**k
            if est < min:
                min = est
        except OverflowError:
            pass
    return min


def HashimotoMQ(n, m, q):
    # Underdetermined case
    minest = 2**1024
    for k in range(0, m):
        for a in range(2, m - k):
            if (n >= ((a + 1) * (m - k - a + 1))) and (n >= (a * (m - k) - (a - 1)**2 + k)):
                h1 = MQ(m - a - k, m - a, q)
                h2 = MQ(a - 1, a - 1, q)
                h3 = MQ(a, a, q)
                est = q**k * (h1 + h2) + (m - a - k + 1) * h3
                if est < minest:
                    minest = est

    return minest


def SolveMQest(n, m, q):
    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    if n > m:
        est = min(HashimotoMQ(n, m, q), HybridMQ(m, m, q))
    else:
        est = HybridMQ(n, m, q)

    return math.floor(math.log2(qcomplex) + math.log2(est))


def reconciliation(v, o, m1, q, l):
    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    if q % 2:
        M = m1 * (l * (l + 1)) // 2
    else:
        M = m1 * l**2

    return math.floor(math.log2(qcomplex) + math.log2(HybridMQ(l * v, M, q)) + logq * max(0, l * v - m1 * l**2))


def wedge(v, o, l, m1, q):
    # using https://gitlab.di.ens.fr/hle/snova_wedge_attack
    # commit 87f695020642701e26c5109666a9c8d06fcecf83
    try:
        from wedge_rank import find_best_projection
    except ModuleNotFoundError:
        print('Import error. Try: export PYTHONPATH=`pwd`', file=sys.stderr)
        quit()

    def cost_complexity(s, v, q, l):
        term1 = 2 * s  # Sparse matrix solver complexity
        term2 = 2 * math.log2(v+1)  # Density (v+1)^2
        L = l * math.log2(q)
        field_cost = 3 * (2 * (L ** 2) + L)
        term3 = math.log2(field_cost)  # Field-arithmetic cost: 3 * (2*(log2(q^l))^2  +  log2(q^l))
        return math.floor(term1 + term2 + term3)

    # Search settings
    STRICT = True          # require S > vars (set False for S >= vars)
    MAX_STATES = 20000     # raise if "not_found_within_budget" occurs often

    # Projection regime:
    #   - "unbalanced" (default):multiset / sorted tuple search
    #   - "balanced": all-equal tuple search
    PROJECTION = "unbalanced"
    n = v + o

    aggr = find_best_projection(
        n, v, m1, l,
        mode="aggressive",
        projection=PROJECTION,
        strict=STRICT,
        max_states=MAX_STATES,
    )

    aggr_est = cost_complexity(math.log2(aggr['vars']), v, q, l) if aggr.get("status") == "found" else 9999

    return aggr_est


def lifting_reconciliation(v, o, q, l, m1):
    try:
        c = math.ceil(v / o)
        if True and q % 2:
            return SolveMQest(c * v, (c * (c + 1)) // 2 * m1, q**l)
        else:
            return SolveMQest(c * v, c**2 * m1, q**l)
    except:
        return 9999


def lifting_intersection(v, o, q, l, m1):
    if q % 2:
        est = HybridMQ(2 * v - o, 3 * m1 - 2, q**l)
    else:
        est = HybridMQ(2 * v - o, 4 * m1 - 2, q**l)

    return math.floor(math.log2(qcomplex) + math.log2(est * q**(l * (v - 2 * o))))


for var in variants:
    try:
        v = var[0]
        o = var[1]
        q = var[2]
        l = var[3]

        if len(var) > 4:
            r = var[4]
            if len(var) > 5:
                m1 = var[5]
            else:
                m1 = math.ceil(o * r / l)
        else:
            r = l
            m1 = math.ceil(o * r / l)

        symmetric = q % 2 != 0

        logq = math.log2(q)
        qcomplex = 2 * logq**2 + logq

        m2 = l * r * o
        n = o + v

        if q == 16:
            PACK_GF = 2
            PACK_BYTES = 1

        elif q == 11:
            PACK_GF = 16
            PACK_BYTES = 7

        elif q == 13:
            PACK_GF = 15
            PACK_BYTES = 7

        elif q == 17:
            PACK_GF = 31
            PACK_BYTES = 16

        elif q == 19:
            PACK_GF = 15
            PACK_BYTES = 8

        elif q == 23:
            PACK_GF = 7
            PACK_BYTES = 4

        elif q == 31:
            PACK_GF = 8
            PACK_BYTES = 5

        elif q == 127:
            PACK_GF = 8
            PACK_BYTES = 7

        elif q == 127**3:
            PACK_GF = 8
            PACK_BYTES = 21

        elif q == 256:
            PACK_GF = 1
            PACK_BYTES = 1

        else:
            raise Exception('Unsupported q')

        def BYTES_GF(x):
            return (PACK_BYTES * (x) + PACK_GF - 1) // PACK_GF

        if symmetric or l == 1:
            pk = BYTES_GF(m1 * o * l * (o * l + 1) // 2) + 16
        else:
            pk = BYTES_GF(m1 * o * l * (o * l)) + 16

        sig = BYTES_GF(r * n * l) + 16

        # Estimates

        cs = math.log2(m1 * n ** 2 * l**2 * r**2 * qcomplex)
        c = (17.5 + cs) / 2 + 1

        logq = math.log2(q)
        qcomplex = 2 * logq**2 + logq

        # Forgery
        drop = l - 1

        direst_est = SolveMQest(l * r * n, m2, q)

        mlc_est = math.floor(math.log2(q) / 2 * m2 + c)

        collision = math.floor(max(logq / 2 * m2 + math.log2(2 * math.sqrt(44 * m2 * 2**17 * 2**35)),
                                   logq * m2 + 15 - 167))

        beullens_est = SolveMQest(l * n - drop, min(m2, m1 * l**2) - drop, q)

        # Key recovery
        kipnisshamir_est = math.floor(math.log2(qcomplex) + logq * l * (v - o))

        reconciliation_est = reconciliation(v, o, m1, q, l)

        if symmetric:
            intersection_est = HybridMQ(l * n + 1, 2 * l**2 * m1 + l * m1 - 2 * l, q)
        else:
            intersection_est = HybridMQ(l * n + 1, 4 * l**2 * m1 - 2 * l, q)

        intersection_est = math.floor(math.log2(qcomplex) + math.log2(intersection_est) +
                                      logq * (l * v - 2 * l * o + 1))

        wedge_est = wedge(v, o, l, m1, q)

        # Lifting
        lifting_ks_est = math.floor(math.log2(2 * logq**2 * l**2 + logq * l) + logq * l * (v - o))

        lifting_rec_est = lifting_reconciliation(v, o, q, l, m1)

        try:
            lifting_int_est = lifting_intersection(v, o, q, l, m1)
        except:
            lifting_int_est = 9999

        # Output
        print(f'(v={v}, o={o}, q={q}, l={l}, r={r}, m1={m1})\t& {pk} & {sig} &', end='')
        print(f'\t\t& {direst_est} & {collision} & {mlc_est} & {beullens_est} &', end='')
        print(f'\t   {reconciliation_est} & MHrec & {kipnisshamir_est} & {intersection_est} & MHint & {wedge_est} &', end='')
        print(f'\t  {lifting_rec_est} & {lifting_ks_est} & {lifting_int_est}', end='')
        print('   \\\\')
        sys.stdout.flush()

    except Exception as exc:
        traceback.print_exc()
