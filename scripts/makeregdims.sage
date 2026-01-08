# Create table of D_reg

import sys

cutoff_degree = 1024
R = PowerSeriesRing(ZZ, 'x', default_prec=cutoff_degree)
x = R.gens()[0]

powzero = 0
for idx in range(cutoff_degree):
    powzero = powzero + x ** idx

print('cutoff_degree =', cutoff_degree)

print('d_reg = {')

for m in range(1, cutoff_degree + 1):
    pol = (1 - x ** 2) ** m
    for n in range(1, m):
        pol = pol / (1 - x)

        num = 1
        for idx in pol.coefficients()[1:]:
            if idx < 1:
                print(' ({}, {}): {},'.format(n, m, num))
                break
            num += 1
    print(' ({}, {}): {},'.format(m, m, m + 1))

    sys.stdout.flush()

print('}')
