# Create table of D_reg

import sys

cutoff_degree = 1024
R = PolynomialRing(ZZ, 'x')
x = R.0

powzero = 0
for idx in range(cutoff_degree):
    powzero = powzero + x ** idx

print('cutoff_degree =', cutoff_degree)

print('d_reg = {')

for m in range(1, cutoff_degree + 1):
    pol = (1 - x ** 2) ** m
    for n in range(1, m + 1):
        pol = pol * powzero

        num = 1
        for idx in pol.coefficients()[1:cutoff_degree]:
            if idx < 1:
                print(' ({}, {}): {},'.format(n, m, num))
                break
            num += 1

    sys.stdout.flush()

print('}')
