# SPDX-License-Identifier: MIT
#
# Script to estimate SNOVA security with respect to Multi-Homogeneous attacks
#
# Copyright (c) 2026 SNOVA TEAM

from itertools import combinations
import traceback
import sys
import math


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


def Reconciliation(v, m1, q, l, lp):
    logq = math.log2(q)
    qcomplex = 2 * l**2 * logq**2 + l * logq

    minest = 9999
    mindim = None
    for k in range(m1 * lp):
        multi_regs = MultiDeg(m1 * lp - k, m1, 2 * m1, q**l, l, lp)

        for dims in multi_regs:
            est = prod([math.comb(m1 * lp - k + dims[i], dims[i]) for i in range(lp)])

            est = 2 * math.log2(est)
            est += math.log2(qcomplex * 3 * (m1 * lp - k + 2)**2)
            est += logq * l * (k + max(v - l * m1, 0))

            if est < minest:
                minest = est
                mindim = dims, k

    return minest, mindim


def mh_reconciliation(v, o, m1, q, l):
    logq = l * math.log2(q)
    qcomplex = 2 * logq**2 + logq

    minest = 9999
    for k in range(0, m1):
        dim = RegDim(m1 - k, m1, q**l)
        if dim < 999:
            est = 3 * (m1 - k + 2)**2 * math.comb(m1 - k + dim, dim)**2
            est = math.log2(qcomplex) + math.log2(est)
            est += logq * (k + max(v - l * m1, 0))

            if est < minest:
                minest = est

    return math.floor(minest)


def Intersection(v, o, m1, q, l, lp):
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


if __name__ == "__main__" or __name__ == "sage.all":

    variants = [
        [24, 5, 19, 4],
        [35, 7, 19, 4],
        [44, 6, 19, 4, 6],

        [37, 17, 19, 2],
        [56, 25, 19, 2],
        [75, 33, 19, 2],

        [20, 4, 31, 4],
        [31, 6, 19, 4],
        [42, 8, 19, 4],
    ]

    try:
        for var in variants:
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

            minrec = 9999
            for lp in range(l, l + 1):
                rec_est, rec_dim = Reconciliation(v, m1, q, l, lp)
                print('Recon.', lp, (v, o, q, l, r, m1), rec_est, rec_dim)
                if rec_est < minrec:
                    minrec = rec_est
                sys.stdout.flush()

            minint = 9999
            for lp in range(l, l + 1):
                int_est, int_dim = Intersection(v, o, m1, q, l, lp)
                print('Int.', lp, '??' if lp * (v + o) > 4 * o * lp**2 - 2 * lp else ' ',
                      (v, o, q, l, r, m1), int_est, int_dim)
                if int_est < minint:
                    minint = int_est
                sys.stdout.flush()

            # mh_rec = mh_reconciliation(v, o, m1, q, l)
            minrec = math.floor(minrec)
            minint = math.floor(minint)

            print()
            print(f'(v={v}, o={o}, q={q}, l={l}, r={r}, m1={m1})\t&  {minrec}   &  {minint}')
            print()
            sys.stdout.flush()

    except Exception as exc:
        traceback.print_exc()
