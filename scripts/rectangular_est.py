# SPDX-License-Identifier: MIT
#
# Script to estimate SNOVA security against various attacks
#
# Copyright (c) 2025 SNOVA TEAM

import math
from regdims import d_reg

params = [
    [24, 5, 23, 4, True],
    [24, 5, 23, 4, False],
    [29, 7, 19, 2, True, 5, 18, 20],
    [29, 7, 19, 2, False, 5, 18, 20],
    [24, 5, 16, 4, True],
    [24, 5, 16, 4, False],
    [43, 17, 16, 2, True],
    [43, 17, 16, 2, False],
    [36, 8, 19, 2, True, 5, 20, 20],
    [36, 8, 19, 2, False, 5, 20, 20],

    [37, 8, 19, 4, True],
    [37, 8, 19, 4, False],
    [52, 8, 19, 2, True, 6, 25, 24],
    [52, 8, 19, 2, False, 6, 25, 24],
    [37, 8, 16, 4, True],
    [37, 8, 16, 4, False],
    [69, 25, 16, 2, True],
    [69, 25, 16, 2, False],

    [60, 10, 23, 4, True],
    [60, 10, 23, 4, False],
    [64, 11, 19, 2, True, 6, 34, 30],
    [64, 11, 19, 2, False, 6, 34, 30],
    [60, 10, 16, 4, True],
    [60, 10, 16, 4, False],
    [99, 33, 16, 2, True],
    [99, 33, 16, 2, False],
]


def GB(n, m):
    D = d_reg[(n, m)]
    return 3 * math.comb(n - 1 + D, D)**2 * math.comb(n + 1, 2)


def HybridMQ(n, m, q):
    # Overdetermined case
    min = GB(m, m)
    for k in range(max(0, n - m), n):
        try:
            est = GB(n - k, m) * q**k
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
                h1 = HybridMQ(m - a - k, m - a, q)
                h2 = HybridMQ(a - 1, a - 1, q)
                h3 = HybridMQ(a, a, q)
                est = q**k * (h1 + h2) + (m - a - k + 1) * h3
                if est < minest:
                    minest = est

    return minest


def SolveMQest(n, m, q):
    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    if n > m:
        est = HashimotoMQ(n, m, q)
    else:
        est = HybridMQ(n, m, q)

    return math.floor(math.log2(qcomplex) + math.log2(est))


def HybridMQest(n, m, q):
    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    return math.floor(math.log2(qcomplex) + math.log2(HybridMQ(n, m, q)))


def lifting_reconciliation(v, o, q, l, m1):
    c = math.ceil(v / o)
    try:
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


for param in params:
    v = param[0]
    o = param[1]
    q = param[2]
    l = param[3]
    aes = param[4]
    if aes:
        continue

    if len(param) > 5:
        r = param[5]
        m1 = param[6]
        Nalpha = param[7]
    else:
        r = l
        m1 = math.ceil(o * r / l)
        Nalpha = r**2 + r

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

    logq = math.log2(q)
    qcomplex = 2 * logq**2 + logq

    # Forgery
    direst_est = SolveMQest(l * r * n, m2, q)

    cs = math.log2(m1 * n ** 2 * l**2 * r**2 * qcomplex)
    c = (17.5 + cs) / 2 + 1
    mlc_est = math.floor(math.log2(q) / 2 * m2 + c)

    # Conservative rank drop
    drop = l + 1
    beullens_est = SolveMQest(l * n - drop, min(m2, m1 * l**2) - drop, q)

    # Key recovery
    kipnisshamir_est = math.floor(math.log2(qcomplex) + logq * l * (v - o))

    if symmetric:
        reconciliation_est = HybridMQest(l * v + 1, m1 * (l * (l + 1)) // 2, q)
    else:
        reconciliation_est = HybridMQest(l * v + 1, m1 * l ** 2, q)

    if symmetric:
        intersection_est = HybridMQ(l * n + 1, 2 * l**2 * m1 + l * m1 - 2 * l, q) * q**(l * v - 2 * l * o + 1)
    else:
        intersection_est = HybridMQ(l * n + 1, 4 * l**2 * m1 - 2 * l, q) * q**(l * v - 2 * l * o + 1)
    intersection_est = math.floor(math.log2(qcomplex) + math.log2(intersection_est))

    # Lifting
    lifting_ks_est = math.floor(math.log2(2 * logq**2 * l**2 + logq * l) + logq * l * (v - o))

    lifting_rec_est = lifting_reconciliation(v, o, q, l, m1)

    lifting_int_est = lifting_intersection(v, o, q, l, m1)

    # Output
    fmt = '{}\t & {} & {} &    \t{} & & {} & {} & \t   {} & & {} & {} & \t  {} & {} & {}'
    print(fmt.format((v, o, q, l, r, m1), pk, sig,
                     direst_est, mlc_est, beullens_est,
                     reconciliation_est, kipnisshamir_est, intersection_est,
                     lifting_rec_est, lifting_ks_est, lifting_int_est))
