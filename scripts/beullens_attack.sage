# SPDX-License-Identifier: MIT
#
# Tool to establish rank of the ABQ construction against Beullens' attack
#
# Copyright (c) 2026 SNOVA TEAM


import math
from hashlib import shake_256
import sys


o = 5
q = 16
l = 4
r = 4

m1 = math.floor(o * r / l)
n_alpha = l * r + 2 * r

m2 = o * l * r


GF_q = GF(q, 'x')
x = GF_q.gen()


if q == 16:
    def from_int(x): return GF_q.from_integer(x)
    def to_int(x): return x.to_integer()
else:
    def from_int(x): return x
    def to_int(x): return int(x)


QR = PolynomialRing(GF_q, 'z')
z = QR.gen()


# Set constants

if q == 11:
    Q_A = 0
    Q_B = 3
    Q_C = 6

    S = matrix(GF_q, l, l)
    for i in range(l):
        for j in range(i, l):
            S[i, j] = (Q_A + i + j) & Q_B
            S[j, i] = S[i, j]
    S[l - 1, l - 1] = Q_C
elif q == 16:
    S = matrix(GF_q, l, l, lambda i, j: from_int(abs(8 - (i + j))))
else:
    raise Exception()

p = S.charpoly('z')
R = companion_matrix(p)

C = [R ^ a for a in range(l)]


def convert_bytes_to_GF(data):
    # Expand public XOF data
    if q == 16:
        res = []
        for item in data:
            res.append(from_int(item % 16))
            res.append(from_int(item // 16))
        return res
    else:
        return [item % q for item in data]


def gen_ABQ(seed):
    # Generate public ABQ

    def create_AB(data, r1, r2):
        # Improve public matrices
        M = matrix(GF_q, r1, r2, lambda i, j: data[i * r2 + j])
        if l == r1 and l == r2:
            f1 = 1
            while M.det() == 0 and f1 < q:
                M += from_int(f1) * S
                f1 += 1
            if f1 == q:
                raise Exception('f1 == q')
        return M

    def nonzero_q(data):
        coefs = [to_int(data[i]) for i in range(l)]
        if coefs[l - 1] == 0:
            coefs[l - 1] = q - (coefs[0] if coefs[0] != 0 else 1)
        return [from_int(coefs[i]) for i in range(l)]

    NUM_ABQ = o * n_alpha * (r * (r + l) + 2 * l)
    abqdata = shake_256(seed).digest(NUM_ABQ)
    abqdata = convert_bytes_to_GF(abqdata)

    A = [create_AB(abqdata[i * r**2:], r, r) for i in range(o * n_alpha)]
    B = [create_AB(abqdata[o * n_alpha * r**2 + i * l * r:], r, l) for i in range(o * n_alpha)]
    q1 = [nonzero_q(abqdata[o * n_alpha * r * (r + l) + i * l:]) for i in range(o * n_alpha)]
    q2 = [nonzero_q(abqdata[o * n_alpha * r * (r + l) + o * n_alpha * l + i * l:]) for i in range(o * n_alpha)]

    return A, B, q1, q2


def get_etilde(r_coeffs):
    etilde = matrix(GF_q, r * l * o, m1 * l**2)

    for mi in range(o):
        for alpha in range(n_alpha):
            mia = mi * n_alpha + alpha
            mi_prime = (mi + alpha) % m1

            # A, q1
            z1a = matrix(GF_q, r, l)
            for i1 in range(r):
                for a1 in range(l):
                    for i2 in range(r):
                        z1a[i1, a1] += A[mia][i1, i2] * r_coeffs[i2, a1]

            z1b = matrix(GF_q, l, l)
            for a1 in range(l):
                for a2 in range(l):
                    for a in range(l):
                        z1b[a1, a2] += C[a1][a2, a] * q1[mia][a]

            z1 = matrix(GF_q, r, l)
            for i1 in range(r):
                for a2 in range(l):
                    for a1 in range(l):
                        z1[i1, a2] += z1a[i1, a1] * z1b[a1, a2]

            # B, q2
            z2a = matrix(GF_q, l, l)
            for j1 in range(l):
                for b1 in range(l):
                    for j2 in range(r):
                        z2a[j1, b1] += B[mia][j2, j1] * r_coeffs[j2, b1]

            z2b = matrix(GF_q, l, l)
            for b1 in range(l):
                for b2 in range(l):
                    for b in range(l):
                        z2b[b1, b2] += C[b1][b2, b] * q2[mia][b]

            z2 = matrix(GF_q, l, l)
            for j1 in range(l):
                for b2 in range(l):
                    for b1 in range(l):
                        z2[j1, b2] += z2a[j1, b1] * z2b[b1, b2]

            for i1 in range(r):
                for j1 in range(l):
                    for a2 in range(l):
                        for b2 in range(l):
                            etilde[(mi * r + i1) * l + j1, (mi_prime * l + a2) * l + b2] += z1[i1, a2] * z2[j1, b2]

    return etilde


index = 0

seed = f'SNOVA_ABQ_{index}'.encode() if index else b'SNOVA_ABQ'

A, B, q1, q2 = gen_ABQ(seed)


for idx in range(0, 100):
    r_c = matrix(GF_q, r, l)
    r_c[0, 0] = 1

    idx1 = idx
    for j1 in range(l):
        for i1 in range(1, r):
            r_c[i1, j1] = from_int(idx1 % q)
            idx1 = idx1 // q

    emat = get_etilde(r_c)
    rank = emat.rank()
    if rank < min(m2, m1 * l**2):
        print('-->', idx, rank, m2, m1 * l**2)
    sys.stdout.flush()
