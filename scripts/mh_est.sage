# SPDX-License-Identifier: MIT
#
# Script to estimate SNOVA security with respect to Multi-Homogeneous attacks
#
# Copyright (c) 2026 SNOVA TEAM

from itertools import combinations
import traceback
import sys
import math


def MultiDeg(np, mp1, mp2, q, l):
    precision = 120 // l
    if l > 1:
        R = PowerSeriesRing(ZZ, l, 't', default_prec=precision)
        t = R.gens()
    else:
        R = PowerSeriesRing(ZZ, 't0', default_prec=precision)
        t = R.gens()

    pol = prod([((1 - t[i] ** 2) / (1 - t[i] ** (2 * q))) ** mp1 *
               ((1 - t[i]**q) / (1 - t[i]))**np / (1 - t[i]) for i in range(l)])
    for i, j in combinations(range(l), 2):
        pol *= (1 - t[i] * t[j]) ** mp2

    multi_regs = {}
    pol_coef = pol.coefficients()
    if l == 1:
        for d_reg in range(precision):
            try:
                if pol_coef[d_reg] <= 0:
                    return [(d_reg,)]
            except:
                pass
        return []

    for item in pol_coef:
        if pol_coef[item] <= 0:
            dsol = reversed(sorted(item.exponents()[0]))
            multi_regs[tuple(dsol)] = 0

    return sorted(multi_regs.keys(), key=lambda x: (sum(x), x))


def Reconciliation(v, m1, q, l, lp):
    logq = math.log2(q)
    qcomplex = 2 * lp**2 * logq**2 + lp * logq

    minest = 9999
    for k in range(0, m1 * l, lp):
        multi_regs = MultiDeg(m1 * l - k, m1, m1 if q % 2 else 2 * m1, q, lp)

        for dims in multi_regs:
            est = prod([math.comb(m1 * l - k + dims[i], dims[i]) for i in range(lp)])

            est = 2 * math.log2(est)
            est += math.log2(qcomplex * 3 * math.comb(m1 * l - k + 2, 2))
            est += logq * (k * lp + l * v - m1 * l**2)

            if est < minest:
                minest = est
                mindim = dims, k

    return minest, mindim


def Intersection(v, o, m1, q, l, lp):
    logq = lp * math.log2(q)
    qcomplex = 2 * logq**2 + logq
    n = o + v

    minest = 9999
    mindim = None

    if q % 2:
        multi_regs = MultiDeg(n, 3 * m1 - 2, 4 * m1, q, lp)
    else:
        multi_regs = MultiDeg(n, 4 * m1 - 2, 8 * m1, q, lp)

    for dims in multi_regs:
        est = prod([math.comb(n + dims[i] - 1, dims[i]) for i in range(lp)])
        est = 2 * math.log2(est)
        est += math.log2(qcomplex * 3 * math.comb(n + 1, 2))
        est += math.log2(q) * l * (v - 2 * o + 1)
        if est < minest:
            minest = est
            mindim = dims

    return minest, mindim


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

try:
    for var in variants:
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

        minrec = 9999
        for lp in range(1, l + 1):
            rec_est, rec_dim = Reconciliation(v, m1, q, l, lp)
            print('Recon.', lp, (v, o, q, l, r, m1), rec_est, rec_dim)
            if rec_est < minrec:
                minrec = rec_est
            sys.stdout.flush()

        minint = 9999
        for lp in range(1, l + 1):
            int_est, int_dim = Intersection(v, o, m1, q, l, lp)
            print('Int.', lp, '??' if lp * (v + o) > 4 * o * lp**2 - 2 * lp else ' ',
                  (v, o, q, l, r, m1), int_est, int_dim)
            if int_est < minint:
                minint = int_est
            sys.stdout.flush()

        print()
        print(tuple(var), math.floor(minrec), math.floor(minint))
        print()
        sys.stdout.flush()

except Exception as exc:
    traceback.print_exc()
