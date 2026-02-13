# SPDX-License-Identifier: MIT
#
# Script to estimate SNOVA security against various attacks
#
# Copyright (c) 2026 SNOVA TEAM


import math
import sys
import traceback

variants = [
    [24, 5, 23, 4],
    [24, 5, 16, 4],
    [38, 7, 19, 2, 5, 16, 20],
    [43, 17, 16, 2],

    # [35, 7, 19, 4],
    [37, 8, 19, 4],
    [37, 8, 16, 4],
    [54, 8, 19, 2, 6, 24, 24],
    [67, 25, 16, 2],

    [60, 10, 23, 4],
    [60, 10, 16, 4],
    [74, 11, 19, 2, 6, 32, 30],
    [90, 33, 16, 2],
]


def RegDim(np, mp, q):
    precision = 200
    R = PowerSeriesRing(ZZ, 't0', default_prec=precision)
    t0 = R.gens()[0]

    pol = ((1 - t0**2) / (1 - t0**(2 * q)))**mp * ((1 - t0**q) / (1 - t0))**np / (1 - t0)
    pol_coef = pol.coefficients()

    for d_reg in range(precision):
        try:
            if pol_coef[d_reg] <= 0:
                return d_reg
        except:
            pass

    return 9999


def GB(n, m, q):
    dim = RegDim(n, m, q)
    est = 3 * math.comb(n + 2, 2) * math.comb(n + dim, dim)**2
    return est


def HybridMQ(n, m, q):
    # Overdetermined case
    min = GB(m, m, q)
    for k in range(max(0, n - m), n):
        try:
            est = GB(n - k, m, q) * q**k
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
            if (n >= ((a + 1) * (m - k - a + 1))) and (n >= (a * (m - k) -
                                                             (a - 1)**2 + k)):
                h1 = GB(m - a - k, m - a, q)
                h2 = GB(a - 1, a - 1, q)
                h3 = GB(a, a, q)
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


def reconciliation(v, m1, q, l):
    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    if q % 2:
        M = m1 * (l * (l + 1)) // 2
    else:
        M = m1 * l**2

    rec1 = math.floor(math.log2(qcomplex) + math.log2(HybridMQ(l * v, M, q)))

    # Using Furue's comment in QR-UOV spec
    rec2 = math.floor(math.log2(qcomplex) + math.log2(HybridMQ(l * v, max(M, l * v), q)))

    return min(rec1, rec2)


def lifting_reconciliation(v, o, q, l, m1):
    try:
        c = math.ceil(v / o)
        if symmetric:
            return SolveMQest(c * v, (c * (c + 1)) // 2 * m1, q**l)
        else:
            return SolveMQest(c * v, c**2 * m1, q**l)
    except:
        return 9999


def lifting_intersection(v, o, q, l, m1):
    if symmetric:
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
            m1 = var[5]
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

        reconciliation_est = reconciliation(v, m1, q, l)

        # Extension field Reconciliation estimate
        ext_reconciliation = math.floor(math.log2(2 * logq**2 * l**2 + logq * l) +
                                        math.log2(HybridMQ(m1 * l, m1 * l, q**l)) + logq*(l * v - l**2 * m1))

        if symmetric:
            intersection_est = HybridMQ(l * n + 1, 2 * l**2 * m1 + l * m1 - 2 * l, q)
        else:
            intersection_est = HybridMQ(l * n + 1, 4 * l**2 * m1 - 2 * l, q)

        intersection_est = math.floor(math.log2(qcomplex) + math.log2(intersection_est) +
                                      logq * (l * v - 2 * l * o + 1))

        # Lifting
        lifting_ks_est = math.floor(math.log2(2 * logq**2 * l**2 + logq * l) + logq * l * (v - o))

        lifting_rec_est = lifting_reconciliation(v, o, q, l, m1)

        try:
            lifting_int_est = lifting_intersection(v, o, q, l, m1)
        except:
            lifting_int_est = 9999

        # Output
        fmt = '{}\t & {} & {} &    \t{} & {} & {} & {} & \t   {} & {} & {} & {} & MHi &\t  {} & {} & {}'
        print(fmt.format((v, o, q, l, r, m1), pk, sig,
                         direst_est, collision, mlc_est, beullens_est,
                         reconciliation_est, ext_reconciliation, kipnisshamir_est, intersection_est,
                         lifting_rec_est, lifting_ks_est, lifting_int_est))
        sys.stdout.flush()

    except Exception as exc:
        traceback.print_exc()
