# SPDX-License-Identifier: MIT
#
# Script to estimate SNOVA security against known attacks
#
# Copyright (c) 2026 SNOVA TEAM


from itertools import combinations
import math
import sys
import traceback

precision = 200
do_mh = True

variants = [
    [26, 5, 19, 4],
    [43, 17, 19, 2],

    [37, 7, 19, 4],
    [64, 25, 19, 2],

    [44, 6, 19, 4, 6],
    [90, 33, 19, 2],
]


plus1 = 0

if plus1:
    print('Using XL')


def regseries(np, mp, pl):
    prec = max(precision, np + 1)
    R = PowerSeriesRing(QQ, 't0', default_prec=prec)
    t0 = R.gens()[0]

    pol = ((1 - t0**2) / (1 - t0**(2 * pl)))**mp * ((1 - t0**pl) / (1 - t0))**np

    return pol.list()


def RegDim(np, mp, q, p1):
    pol_coef = regseries(np + p1, mp, q)
    for d_reg in range(len(pol_coef)):
        if pol_coef[d_reg] <= 0:
            return d_reg
    return math.inf


def MQ(n, m, q, p1):
    dim = RegDim(n, m, q, p1)
    if dim < 999:
        est = 3 * math.comb(n + p1 + 1, 2) * math.comb(n + p1 - 1 + dim, dim)**2
        return est
    else:
        return math.inf


def HybridMQ(n, m, q, p1=plus1):
    # Overdetermined case, n < m
    if m < n:
        raise Exception('HybridMQ', n, m)

    minest = math.inf
    for k in range(n):
        try:
            est = MQ(n - k, m, q, p1) * q**k
            if est < minest:
                minest = est
        except OverflowError:
            pass
    return minest


def HashimotoMQ(n, m, q, p1):
    # Underdetermined case, n > m
    if n < m:
        raise Exception('HashimotoMQ', n, m)

    minest = math.inf
    for k in range(0, m):
        for a in range(2, m - k):
            try:
                if (n >= ((a + 1) * (m - k - a + 1))) and (n >= (a * (m - k) - (a - 1)**2 + k)):
                    h1 = MQ(m - a - k, m - a, q, p1)
                    h2 = MQ(a - 1, a - 1, q, p1)
                    h3 = MQ(a, a, q, p1)
                    est = q**k * (h1 + h2) + (m - a - k + 1) * h3
                    if est < minest:
                        minest = est
            except:
                pass

    return minest


def SolveMQest(n, m, q, p1=plus1):
    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    if n > m:
        est = min(HashimotoMQ(n, m, q, p1), HybridMQ(m, m, q, p1))
    else:
        est = HybridMQ(n, m, q, p1)

    return math.log2(qcomplex) + math.log2(est)


def reconciliation(v, o, m1, q, l):
    logq = math.log2(q)

    if q % 2:
        M = m1 * (l * (l + 1)) // 2
    else:
        M = m1 * l**2

    return SolveMQest(l * v, M, q) + logq * max(0, l * v - m1 * l**2)


def pl_T(aa, pl, d):
    sum = 0

    for i in range(min(aa, math.floor(d / pl)) + 1):
        sum += (-1)**i * math.comb(aa, i) * math.comb(d - i * pl + aa - 1, aa - 1)

    return sum


def pl_reconciliation(p, e, N, M, a):
    minest = math.inf
    minparam = None

    for l in range(e + 1):
        pl = p**l

        for k in range(a):

            dmin = 0
            nn = N - k
            aa = a - k

            coeffs = regseries(nn, M, pl)

            for d in range(2, aa * (pl - 1) + 1):
                tval = pl_T(aa, pl, d)

                if d >= len(coeffs) or coeffs[d] <= tval:
                    dmin = d
                    break

            if dmin > 0 and dmin < precision:
                est = float(3 * math.comb(nn + 1, 2) * (pl_T(nn, pl, dmin) - pl_T(aa, pl, dmin) + 1)**2)
                if est == 0:
                    continue

                logq = math.log2(p**e)
                qcomplex = 2 * logq**2 + logq
                if est < math.inf:
                    est = math.log2(qcomplex) + math.log2(est)

                if est < minest:
                    minest = est
                    minparam = l, k, dmin

            if dmin >= precision:
                pass
                # print('\tOverflow', l, k, est, dmin)

    return minest


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
        return term1 + term2 + term3

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

    aggr_est = cost_complexity(math.log2(aggr['vars']), v, q, l) if aggr.get("status") == "found" else math.inf

    return aggr_est


def oddwedge(v, o, l, m1, q):
    # using https://gitlab.di.ens.fr/hle/snova_wedge_attack
    # commit 87f695020642701e26c5109666a9c8d06fcecf83
    try:
        from wedge_odd_rank import find_best_projection_odd
    except ModuleNotFoundError:
        print('Import error. Try: export PYTHONPATH=`pwd`', file=sys.stderr)
        quit()

    def cost_complexity(s, v, q, l):
        term1 = 2 * s  # Sparse matrix solver complexity
        term2 = 2 * math.log2(v+1)  # Density (v+1)^2
        L = l * math.log2(q)
        field_cost = 3 * (2 * (L ** 2) + L)
        term3 = math.log2(field_cost)  # Field-arithmetic cost: 3 * (2*(log2(q^l))^2  +  log2(q^l))
        return term1 + term2 + term3

    # Search settings
    STRICT = True          # require S > vars (set False for S >= vars)
    MAX_STATES = 20000     # raise if "not_found_within_budget" occurs often

    # Projection regime:
    #   - "unbalanced" (default):multiset / sorted tuple search
    #   - "balanced": all-equal tuple search
    PROJECTION = "unbalanced"
    n = v + o

    aggr = find_best_projection_odd(
        n, v, m1, l,
        mode="conservative",
        projection=PROJECTION,
        strict=STRICT,
        max_states=MAX_STATES,
    )

    aggr_est = cost_complexity(math.log2(aggr['vars']), v, q, l) if aggr.get("status") == "found" else math.inf

    return aggr_est


# Lars' attack
def ranest(v, o, l, m1, q):
    v1 = v * l
    o1 = o * l

    if q % 2:
        m = m1 * (l * (l - 1)) // 2
    else:
        m = m1 * l**2

    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    def kerndim(v, o, m):
        for op in range(math.floor(2 * v / m + 1), o + 1):
            sum = 0
            sign = 1
            for i in range(0, math.floor(op / 2) + 1):
                sum += sign * math.comb(m + i - 1, i) * math.comb(v + op, v + 2 * i)
                sign = -sign
            # print(op, sum)
            if sum <= 1:
                return op

        return 0

    op = kerndim(v1, o1, m)
    if op == 0:
        return math.inf

    est = 3 * math.comb(v1 + 2, 2) * math.comb(v1 + op, v1)**2

    return math.log2(qcomplex) + math.log2(est)


def lifting_reconciliation(v, o, q, l, m1):
    c = math.ceil(v / o)
    if q % 2:
        # Let's be conservative
        est1 = SolveMQest(c * v, (c * (c + 1)) // 2 * m1, q**l)
        est2 = SolveMQest(c * v, c**2 * m1, q**l)
        return min(est1, est2)
    else:
        return SolveMQest(c * v, c**2 * m1, q**l)


def lifting_intersection(v, o, q, l, m1):

    if v < 2 * o:
        est = SolveMQest(2 * v - o, 3 * m1 - 2, q**l)
        return est

    if q % 2:
        est = SolveMQest(2 * v - o, 3 * m1 - 2, q**l)
    else:
        est = SolveMQest(2 * v - o, 4 * m1 - 2, q**l)
    est = est + math.log2(q) * (l * (v - 2 * o))

    return est


def MultiDeg(np, mp1, mp2, q, l, lp):
    precision = max(120 // l, np * l // lp + 1)
    if lp > 1:
        R = PowerSeriesRing(ZZ, lp, 't', default_prec=precision)
        t = R.gens()
    else:
        R = PowerSeriesRing(ZZ, 't0', default_prec=precision)
        t = R.gens()

    if False:
        pol = prod([((1 - t[i]**2) / (1 - t[i]**(2 * q))) ** mp1 *
                    ((1 - t[i]**q) / (1 - t[i]))**np / (1 - t[i]) for i in range(lp)])
    else:
        pol = prod([(1 - t[i]**2) ** mp1 / (1 - t[i])**(np + 1) for i in range(lp)])

    for i, j in combinations(range(lp), 2):
        pol *= (1 - t[i] * t[j]) ** mp2

    if lp == 1:
        pol_coef = pol.list()
        for d_reg in range(len(pol_coef)):
            if pol_coef[d_reg] <= 0:
                return [(d_reg,)]
        return []

    multi_regs = {}
    pol_coef = pol.coefficients()
    for item in pol_coef:
        if pol_coef[item] <= 0:
            dsol = reversed(sorted(item.exponents()[0]))
            multi_regs[tuple(dsol)] = 0

    return sorted(multi_regs.keys(), key=lambda x: (sum(x), x))


def mh_reconciliation(v, m1, q, l, lp):
    logq = math.log2(q)
    qcomplex = 2 * l**2 * logq**2 + l * logq

    minest = 9999
    mindim = None
    for k in range(m1 * lp):
        try:
            multi_regs = MultiDeg(m1 * lp - k, m1, 2 * m1, q**l, l, lp)
        except:
            print('MultiDeg')
            continue

        for dims in multi_regs:
            est = prod([math.comb(m1 * lp - k + dims[i], dims[i]) for i in range(lp)])

            est = 2 * math.log2(est)
            est += math.log2(qcomplex * 3 * (m1 * lp - k + 2)**2)
            est += logq * l * (k + max(v - l * m1, 0))

            if est < minest:
                minest = est
                mindim = dims, k

    return minest, mindim


def mh_intersection(v, o, m1, q, l, lp):
    logq = l * math.log2(q)
    qcomplex = 2 * logq**2 + logq
    n = o + v

    minest = 9999
    mindim = None

    if q % 2:
        multi_regs = MultiDeg(n, 3 * m1 - 2, 4 * m1, q, l, lp)
    else:
        multi_regs = MultiDeg(n, 4 * m1 - 2, 8 * m1, q, l, lp)

    for dims in multi_regs:
        est = prod([math.comb(n + dims[i] - 1, dims[i]) for i in range(lp)])
        est = 2 * math.log2(est)
        est += math.log2(qcomplex * 3 * math.comb(n + 1, 2))
        est += math.log2(q) * l * (v - 2 * o + 1)
        if est < minest:
            minest = est
            mindim = dims

    return minest, mindim


for var in variants:
    if len(var) == 0:
        print()
        continue

    try:
        v = var[0]
        o = var[1]
        q = var[2]
        l = var[3]
        m2 = 0

        if len(var) > 4:
            r = var[4]
            if len(var) > 5:
                m1 = var[5]
                if len(var) > 6:
                    m2 = var[6]
            else:
                m1 = math.ceil(o * r / l)
        else:
            r = l
            m1 = math.ceil(o * r / l)

        symmetric = q % 2 != 0

        logq = math.log2(q)
        qcomplex = 2 * logq**2 + logq

        if m2 == 0:
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

        elif q == 32:
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

        cs = math.log2(m1 * n**2 * l**3 * r * qcomplex)
        c = (17.5 + cs) / 2 + 1

        logq = math.log2(q)
        qcomplex = 2 * logq**2 + logq

        # Forgery

        direct_est = SolveMQest(l * r * n, m2, q, p1=1)

        mlc_est = math.log2(q) / 2 * m2 + c

        collision = max(logq / 2 * m2 + math.log2(2 * math.sqrt(44 * m2 * 2**17 * 2**35)),
                        logq * m2 + 15 - 167)

        drop = l - 1
        beullens_est = SolveMQest(l * n - drop, min(m2, m1 * l**2) - drop, q, p1=1)

        # Key recovery

        if len(sys.argv) == 7 or l == 1:
            kipnisshamir_est = math.log2(qcomplex) + logq * l * (v - o) + 2.8 * math.log2(n)
        else:
            # Lower bounded by
            kipnisshamir_est = math.log2(qcomplex) + logq * l * (v - o)

        reconciliation_est = reconciliation(v, o, m1, q, l)

        if q == 16:
            pl_reconciliation_est = pl_reconciliation(2, 4, (o + v) * l, m1 * l**2, o * l)
        else:
            pl_reconciliation_est = pl_reconciliation(q, 1, (o + v) * l, m1 * l**2, o * l)

        if v < 2 * o:
            if symmetric:
                intersection_est = SolveMQest(l * (2 * v - o) + 1, 2 * l**2 * m1 + l * m1 - 2 * l, q)
            else:
                intersection_est = SolveMQest(l * (2 * v - o) + 1, 4 * l**2 * m1 - 2 * l, q)
        else:
            if symmetric:
                intersection_est = SolveMQest(l * n + 1, 2 * l**2 * m1 + l * m1 - 2 * l, q)
            else:
                intersection_est = SolveMQest(l * n + 1, 4 * l**2 * m1 - 2 * l, q)
            intersection_est = intersection_est + logq * (l * v - 2 * l * o + 1)

        if q % 2:
            wedge_est = oddwedge(v, o, l, m1, q)
        else:
            wedge_est = wedge(v, o, l, m1, q)

        rans = ranest(v, o, l, m1, q)

        # Lifting
        lifting_ks_est = math.log2(2 * logq**2 * l**2 + logq * l) + logq * l * (v - o)

        lifting_rec_est = lifting_reconciliation(v, o, q, l, m1)

        lifting_int_est = lifting_intersection(v, o, q, l, m1)

        # Multi Homogeneous
        if do_mh:
            mh_rec, _ = mh_reconciliation(v, m1, q, l, l)
            mh_int, _ = mh_intersection(v, o, m1, q, l, l)
        else:
            mh_rec = math.inf
            mh_int = math.inf

        # Output
        def f(x): return math.inf if x == math.inf else math.floor(x)

        direct_est = f(direct_est)
        collision = f(collision)
        mlc_est = f(mlc_est)
        beullens_est = f(beullens_est)
        reconciliation_est = f(reconciliation_est)
        mh_rec = f(mh_rec)
        kipnisshamir_est = f(kipnisshamir_est)
        intersection_est = f(intersection_est)
        mh_int = f(mh_int)
        wedge_est = f(wedge_est)
        pl_reconciliation_est = f(pl_reconciliation_est)
        rans = f(rans)
        lifting_rec_est = f(lifting_rec_est)
        lifting_ks_est = f(lifting_ks_est)
        lifting_int_est = f(lifting_int_est)

        print(f'(v={v}, o={o}, q={q}, l={l}, r={r}, m1={m1})\t& {pk} & {sig} &', end='')
        print(f'\t\t& {direct_est} & {collision} & {mlc_est} & {beullens_est} &', end='')
        print(f'\t   {reconciliation_est} & {mh_rec} & {kipnisshamir_est} & {intersection_est} & {mh_int} & {wedge_est} & {pl_reconciliation_est} &', end='')
        print(f'\t  {lifting_rec_est} & {lifting_ks_est} & {lifting_int_est}', end='')
        print('   \\\\')
        sys.stdout.flush()

    except Exception as exc:
        traceback.print_exc()
