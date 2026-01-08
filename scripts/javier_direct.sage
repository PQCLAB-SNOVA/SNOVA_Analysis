'''
Analysis script by Javier Verbel (October 30, 2024)
https://eprint.iacr.org/2024/1770

Edited by SNOVA team 2026
'''

from sage.all_cmdline import *

from itertools import combinations, product
from math import log2, floor, comb as binomial
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ

import math
import sys
import datetime


def multi_degree(mon):
    return mon.exponents()[0]


def generate_bounded_combinations(multi_deg):
    # For each d_i in the vector d, generate a range from 0 to d_i
    ranges = [range(di + 1) for di in multi_deg]

    # Use itertools.product to generate the Cartesian product of these ranges
    combinations = list(product(*ranges))
    return combinations[:-1]  # movind multi_deg that should be in the last column


def n_mon(multi_deg, n_vars_per_set):
    output = 1
    for a in multi_deg:
        output *= binomial(n_vars_per_set + a - 1, a)
    return output


def num_columns(n_vars_per_set, p, multi_deg):
    ncols_ = 0
    smaller_multi_degrees = generate_bounded_combinations(multi_deg)
    for vec_d in smaller_multi_degrees:
        ncols_ += n_mon(vec_d, n_vars_per_set)

    ncols_ *= p
    ncols_ += n_mon(multi_deg, n_vars_per_set)
    return ncols_


def special_xl_complexity(n_vars_per_set, p, multi_deg):
    if n_vars_per_set == 0:
        time = 0
    else:
        time = log2(3 * n_vars_per_set ** 2) + 2 * log2(num_columns(n_vars_per_set, p, multi_deg))
    return time


def lry_series(o, l, k=0, prec=None):
    if l == 2:
        names = 't1,t2'
    elif l == 3:
        names = 't1,t2,t3'
    elif l == 4:
        names = 't1,t2,t3,t4'
    else:
        raise ValueError

    if prec is None:
        prec = o * l * r + 2
        if l == 3:
            prec = 20
        elif l == 4:
            prec = 20
    # print(prec)
    PS = PowerSeriesRing(QQ, names, default_prec=prec)
    t = PS.gens()
    n_i = l * [r * o - k // l]
    series = 1
    for i, j in combinations(range(l), 2):
        series *= (1 - t[i] * t[j]) ** (2 * o)
    for i in range(l):
        series *= (1 - t[i] ** 2) ** (o)
    for i in range(l):
        series /= (1 - t[i]) ** (n_i[i])
    return series


def ffd_estimation(o, l, p, k=0, prec=None):
    series = lry_series(o, l, k=k, prec=prec)
    coef_dic = series.coefficients()
    t = series.parent().gens()
    monomials_in_series = coef_dic.keys()
    for mon, coef in coef_dic.items():
        if coef < 0:
            sum_ = coef
            target_multi_degree = multi_degree(mon)
            smaller_multi_degrees = generate_bounded_combinations(target_multi_degree)
            for current_multi_degree in smaller_multi_degrees:
                current_mon = prod([t[i] ** current_multi_degree[i] for i in range(l)])
                if current_mon in monomials_in_series:
                    sum_ += p * coef_dic[current_mon]
            if sum_ < 0:
                # print(target_multi_degree, sum_)
                return target_multi_degree


def complexity_find_z(v, o, q, l, rank_tilde_E):
    return (o * l * r - rank_tilde_E) * log2(q) + 6 * log2(l)


def complexity_given_attack_params(v, o, q, l, rank_tilde_E, value_k, multi_deg):
    p = o * l * r - rank_tilde_E
    n_vars_per_set = o * r - value_k // l
    time_t = (value_k - p) * log2(q) + special_xl_complexity(n_vars_per_set, p, multi_deg)
    time_t += log2(2 * log2(q ** l) ** 2 + log2(q ** l))  # Adding bit-complexity factor
    return time_t


def new_direct_complexity(v, o, q, l, rank_tilde_E, range_k=None):
    # print(v, o, q, l)
    sys.stdout.flush()
    # "Not the best way to find the optimal. TOBEUPDATED"
    p = o * l * r - rank_tilde_E
    optimal_k = 0
    opt_multi_deg = 0
    time = 1000
    time_z = complexity_find_z(v=v, o=o, q=q, l=l, rank_tilde_E=rank_tilde_E)
    if range_k is None:
        range_k = range(l, l * o * r, l)
    for value_k in range_k:
        if value_k > p:
            multi_deg = ffd_estimation(o, l, p, k=value_k)
            if multi_deg is None:
                continue
            t_time = complexity_given_attack_params(v, o, q, l, rank_tilde_E, value_k, multi_deg)
            # print(datetime.datetime.now(), value_k, multi_deg, t_time, time_z)
            sys.stdout.flush()
            if t_time < time:
                time = t_time
                optimal_k = value_k
                opt_multi_deg = multi_deg

    return max(time_z, time), time_z, time, optimal_k, opt_multi_deg


variants = [
    # SL 1 RectSNOVA
    [0, 29, 7, 2, 5, 19],
    [0, 36, 8, 2, 5, 19],

    # SL 3 RectSNOVA
    [0, 46, 9, 2, 5, 23],
    [0, 52, 8, 2, 6, 19, 25],

    # SL 5 RectSNOVA
    [0, 64, 11, 2, 6, 19, 34],
]

variants += [
    [0, 24, 5, 4, 4, 23],
    [0, 24, 5, 4, 4, 16],
    [0, 43, 17, 2, 2, 16],

    [0, 37, 8, 4, 4, 19],
    [0, 37, 8, 4, 4, 16],
    [0, 69, 25, 2, 2, 16],

    [0, 99, 33, 2, 2, 16],
    [0, 60, 10, 4, 4, 23],
    [0, 60, 10, 4, 4, 16],
]

for var in variants:
    v = var[1]
    o = var[2]
    l = var[3]
    r = var[4]
    q = var[5]

    if len(var) > 6:
        m1 = var[6]
    else:
        m1 = math.ceil(o * r / l)

    for delta in range(10):
        if (m1 * l**2 - o * r * l + delta) * delta > l * (r - 1):
            break
    delta = delta - 1

    if q == 16 and l > 3:
        delta = l + 1

    p = delta
    rank_tilde_E = o * l * r - p
    output = new_direct_complexity(v, o, q, l, rank_tilde_E, range_k=range(l, o * l * r, l))

    print(var, [m1, delta], rank_tilde_E, output, datetime.datetime.now(),
          int(output[0] + l * abs(v - o * r) * log2(q)))
    sys.stdout.flush()
