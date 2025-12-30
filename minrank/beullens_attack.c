/**
 * Analyze rank and minrank of E matrix
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fips202.h"

#if 1

#define SNOVA_o 7
#define SNOVA_q 19
#define SNOVA_l 2
#define SNOVA_r 5
#define SNOVA_m1 18
#define SNOVA_alpha 20

#elif 0

#define SNOVA_o 8
#define SNOVA_q 19
#define SNOVA_l 2
#define SNOVA_r 5
#define SNOVA_m1 20
#define SNOVA_alpha 20

#elif 0

#define SNOVA_o 9
#define SNOVA_q 23
#define SNOVA_l 2
#define SNOVA_r 5
#define SNOVA_m1 23
#define SNOVA_alpha 22

#elif 0

#define SNOVA_o 8
#define SNOVA_q 19
#define SNOVA_l 2
#define SNOVA_r 6
#define SNOVA_m1 25
#define SNOVA_alpha 24

#elif 0

#define SNOVA_o 11
#define SNOVA_q 19
#define SNOVA_l 2
#define SNOVA_r 6
#define SNOVA_m1 34
#define SNOVA_alpha 30

#elif 0

#define SNOVA_o 6
#define SNOVA_q 19
#define SNOVA_l 3
#define SNOVA_r 4
#define SNOVA_m1 9
#define SNOVA_alpha 20

#elif 0

#define SNOVA_o 5
#define SNOVA_q 23
#define SNOVA_l 4
#define SNOVA_r 4
#define SNOVA_m1 5
#define SNOVA_alpha 20

#else

#define SNOVA_o 8
#define SNOVA_q 19
#define SNOVA_l 2
#define SNOVA_r 4
#define SNOVA_m1 17
#define SNOVA_alpha 16

#endif

#ifndef SNOVA_m1
#define SNOVA_m1 ((SNOVA_o * SNOVA_r + SNOVA_l - 1) / SNOVA_l)
#endif

#ifndef SNOVA_alpha
#define SNOVA_alpha (SNOVA_r * SNOVA_r + SNOVA_r)
#endif

#define SNOVA_l2 (SNOVA_l * SNOVA_l)
#define SNOVA_r2 (SNOVA_r * SNOVA_r)
#define SNOVA_lr (SNOVA_l * SNOVA_r)
#define SNOVA_olr16 ((SNOVA_o * SNOVA_lr + 15) / 16)
#define SNOVA_olr (SNOVA_olr16 * 16)

#define SNOVA_m1l2 (((SNOVA_m1 * SNOVA_l2 + 15) / 16) * 16)

#define EMAT_COLS (SNOVA_m1 * SNOVA_r * SNOVA_l * SNOVA_r * SNOVA_l)

uint16_t invtab[SNOVA_q] = {0};
uint16_t mininvtab[SNOVA_q] = {0};

/**
 * Powers of companion matrix to characteristic polynomial of S.
 */
#if SNOVA_l == 2 && SNOVA_q == 19
uint16_t C[2][2][2] = {{{1, 0}, {0, 1}}, {{0, 8}, {1, 16}}};
#elif SNOVA_l == 3 && SNOVA_q == 19
uint16_t C[3][3][3] = {
	{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{0, 0, 15}, {1, 0, 7}, {0, 1, 0}}, {{0, 15, 0}, {0, 7, 15}, {1, 0, 7}}
};
#elif SNOVA_l == 4 && SNOVA_q == 19
uint16_t C[4][4][4] = {{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}},
	{{0, 0, 0, 14}, {1, 0, 0, 16}, {0, 1, 0, 12}, {0, 0, 1, 1}},
	{{0, 0, 14, 14}, {0, 0, 16, 11}, {1, 0, 12, 9}, {0, 1, 1, 13}},
	{{0, 14, 14, 11}, {0, 16, 11, 13}, {0, 12, 9, 15}, {1, 1, 13, 3}}
};
#elif SNOVA_l == 2 && SNOVA_q == 23
uint16_t C[2][2][2] = {{{1, 0}, {0, 1}}, {{0, 5}, {1, 0}}};
#elif SNOVA_l == 3 && SNOVA_q == 23
uint16_t C[3][3][3] = {
	{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{0, 0, 20}, {1, 0, 14}, {0, 1, 3}}, {{0, 20, 14}, {0, 14, 16}, {1, 3, 0}}
};
#elif SNOVA_l == 4 && SNOVA_q == 23
uint16_t C[4][4][4] = {{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}},
	{{0, 0, 0, 7}, {1, 0, 0, 6}, {0, 1, 0, 16}, {0, 0, 1, 4}},
	{{0, 0, 7, 5}, {0, 0, 6, 8}, {1, 0, 16, 1}, {0, 1, 4, 9}},
	{{0, 7, 5, 17}, {0, 6, 8, 13}, {0, 16, 1, 14}, {1, 4, 9, 14}}
};
#else
#error "Unsupported parameters"
#endif

uint64_t get_nsec(void) {
#if 1
	struct timespec time = {0};
	timespec_get(&time, TIME_UTC);
	return (uint64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#else
	// Cycles
	uint32_t lo = 0, hi = 0;
	uint64_t o = 0;
	__asm__ __volatile__("rdtscp" : "=a"(lo), "=d"(hi) : : "%ecx");
	o = hi;
	o <<= 32;
	return (o | lo);
#endif
}

int get_rank(uint16_t matrix[SNOVA_m1 * SNOVA_l2][SNOVA_olr]) {
	int rank = 0;
	for (int j1 = 0; j1 < SNOVA_o * SNOVA_lr; j1++) {
		int i1;

		for (i1 = rank; i1 < SNOVA_m1 * SNOVA_l2 - 1; i1++) {
			if (matrix[i1][j1]) {
				break;
			}
		}

		if (matrix[i1][j1] == 0) {
			continue;
		}

		if (i1 > rank) {
			for (int j2 = 0; j2 < SNOVA_olr; j2++) {
				uint16_t temp = matrix[rank][j2];
				matrix[rank][j2] = matrix[i1][j2];
				matrix[i1][j2] = temp;
			}
		}

		uint16_t minverse = mininvtab[matrix[rank][j1]];

		for (int i2 = rank + 1; i2 < SNOVA_m1 * SNOVA_l2; i2++) {
			uint16_t gji = (minverse * matrix[i2][j1]) % SNOVA_q;
			for (int j2 = 0; j2 < SNOVA_olr; j2++) {
				matrix[i2][j2] += gji * matrix[rank][j2];
			}
		}

		for (int i2 = rank + 1; i2 < SNOVA_m1 * SNOVA_l2; i2++) {
			for (int j2 = 0; j2 < SNOVA_olr; j2++) {
				matrix[i2][j2] = matrix[i2][j2] % SNOVA_q;
			}
		}

		rank++;
	}

	return rank;
}

void get_etilde(uint16_t etilde[SNOVA_m1l2][SNOVA_olr], uint16_t rcoefs[SNOVA_lr], uint16_t emat[EMAT_COLS][SNOVA_olr]) {
	// Helper matrix M
	uint16_t mr[SNOVA_l][SNOVA_lr] = {0};
	for (int i1 = 0; i1 < SNOVA_r; i1++)
		for (int a = 0; a < SNOVA_l; a++)
			for (int a2 = 0; a2 < SNOVA_l; a2++) {
				for (int a1 = 0; a1 < SNOVA_l; a1++) {
					mr[a2][i1 * SNOVA_l + a] += rcoefs[i1 * SNOVA_l + a1] * C[a1][a2][a];
				}
				mr[a2][i1 * SNOVA_l + a] = mr[a2][i1 * SNOVA_l + a] % SNOVA_q;
			}

	// Get E tilde
	for (int ip = 0; ip < SNOVA_m1; ip++)
		for (int a2 = 0; a2 < SNOVA_l; a2++)
			for (int b2 = 0; b2 < SNOVA_l; b2++) {
				uint16_t sum_a[SNOVA_olr] = {0};

				for (int a = 0; a < SNOVA_l; a++)
					for (int i2 = 0; i2 < SNOVA_r; i2++) {
						uint16_t sum_b[SNOVA_olr] = {0};

						for (int j2 = 0; j2 < SNOVA_r; j2++)
							for (int b = 0; b < SNOVA_l; b++)
								for (int k1 = 0; k1 < SNOVA_olr; k1++)
									sum_b[k1] += mr[b2][j2 * SNOVA_l + b] *
									             emat[(((ip * SNOVA_r + i2) * SNOVA_l + a) * SNOVA_r + j2) * SNOVA_l + b][k1];

						for (int k1 = 0; k1 < SNOVA_olr; k1++) {
							sum_b[k1] = sum_b[k1] % SNOVA_q;
						}

						for (int k1 = 0; k1 < SNOVA_olr; k1++) {
							sum_a[k1] += mr[a2][i2 * SNOVA_l + a] * sum_b[k1];
						}
					}

				for (int k1 = 0; k1 < SNOVA_olr; k1++) {
					etilde[(ip * SNOVA_l + a2) * SNOVA_l + b2][k1] = sum_a[k1] % SNOVA_q;
				}
			}
}

// Random permutation
void make_idx_table(uint32_t* idxt, char* seed) {
	uint32_t prng_output_public[EMAT_COLS];

	shake256((uint8_t*)prng_output_public, sizeof(prng_output_public), (uint8_t*)seed, strlen(seed));

	for (int mi = 0; mi < EMAT_COLS; ++mi) {
		idxt[mi] = mi;
	}

	for (int mi = 0; mi < EMAT_COLS - 1; ++mi) {
		int mj = mi + (prng_output_public[mi] % (EMAT_COLS - mi));
		int swap = idxt[mi];
		idxt[mi] = idxt[mj];
		idxt[mj] = swap;
	}
}

#if SNOVA_q == 19
#define Q_A 1
#define Q_B 3
#define Q_C 15
#elif SNOVA_q == 23
#define Q_A 1
#define Q_B 11
#define Q_C 22
#else
#error
#endif

uint16_t gf_S[SNOVA_l * SNOVA_l2] = {0};

static void gen_S_array(void) {
	memset(gf_S, 0, sizeof(gf_S));

	for (int i1 = 0; i1 < SNOVA_l; i1++) {
		gf_S[i1 * SNOVA_l + i1] = 1;
	}

#if SNOVA_l > 1
	// Set S^1, the irreducible S matrix
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_S[SNOVA_l2 + i1 * SNOVA_l + j1] = ((Q_A + i1 + j1) & Q_B) % SNOVA_q;
		}
	gf_S[2 * SNOVA_l2 - 1] = Q_C % SNOVA_q;

	for (int si = 2; si < SNOVA_l; si++) {
		for (int i1 = 0; i1 < SNOVA_l; i1++)
			for (int j1 = 0; j1 < SNOVA_l; j1++) {
				uint16_t sum = 0;
				for (int k1 = 0; k1 < SNOVA_l; k1++) {
					sum += gf_S[SNOVA_l2 + i1 * SNOVA_l + k1] * gf_S[(si - 1) * SNOVA_l2 + k1 * SNOVA_l + j1] % SNOVA_q;
				}
				gf_S[si * SNOVA_l2 + i1 * SNOVA_l + j1] = sum % SNOVA_q;
			}
	}
#endif
}

static uint16_t gf_mat_det(uint16_t* a) {
#define DET_SUB(a, b) (a - b)
#define DET_MULT(a, b) (a * b)
	int32_t det = 0;
#if SNOVA_l == 1
	det = a[0];
#elif SNOVA_l == 2
	det = DET_SUB(DET_MULT(a[0], a[3]), DET_MULT(a[1], a[2]));
#elif SNOVA_l == 3
	det = DET_MULT(a[0], DET_SUB(DET_MULT(a[4], a[8]), DET_MULT(a[5], a[7])));
	det += DET_MULT(a[1], DET_SUB(DET_MULT(a[5], a[6]), DET_MULT(a[3], a[8])));
	det += DET_MULT(a[2], DET_SUB(DET_MULT(a[3], a[7]), DET_MULT(a[4], a[6])));
#elif SNOVA_l == 4
	int32_t DET_l;
	int32_t DET_r;
#define DET_L(x, y) DET_l = DET_SUB(DET_MULT(a[x], a[4 + y]), DET_MULT(a[y], a[4 + x]))
#define DET_R(x, y) DET_r = DET_SUB(DET_MULT(a[8 + x], a[12 + y]), DET_MULT(a[8 + y], a[12 + x]))
#define DET22(x1, y1, x2, y2) \
    DET_L(x1, y1);            \
    DET_R(x2, y2);            \
    det += DET_MULT(DET_l, DET_r)
	DET22(0, 1, 2, 3);
	DET22(0, 2, 3, 1);
	DET22(0, 3, 1, 2);
	DET22(1, 2, 0, 3);
	DET22(1, 3, 2, 0);
	DET22(2, 3, 0, 1);
#undef DET_R
#undef DET22
#undef DET_L
#endif
	return det % SNOVA_q;
}

/**
 * Ensure that a matrix is invertible by adding multiples of S
 */
static inline void be_invertible_by_add_aS(uint16_t* mat, const uint16_t* orig, const int l1, const int l2) {
	memcpy(mat, orig, 2 * l1 * l2);
	if ((l1 == SNOVA_l) && (l2 == SNOVA_l))
		if (gf_mat_det(mat) == 0) {
			for (uint16_t f1 = 1; f1 < SNOVA_q; f1++) {
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					mat[i1] = (mat[i1] + (f1 * gf_S[SNOVA_l2 + i1])) % SNOVA_q;
				}
				if (gf_mat_det(mat) != 0) {
					break;
				}
			}
		}
}

int main(int argc, char** argv) {
	if (argc < 3 || argc > 4) {
		printf("Supply index, partitions, and (optional) number of loops\n");
		return 0;
	}

	gen_S_array();

	uint64_t num_r = 1;
	uint64_t start = 0;
	uint64_t number_of_keys;

	int64_t num = 0;
	int64_t tot = 1;
	int64_t loops = 1;

	sscanf(argv[1], "%li", &num);
	sscanf(argv[2], "%li", &tot);
	if (num > tot) {
		printf("Error: index > number of partitions\n");
		return 0;
	}
	if (argc == 4) {
		sscanf(argv[3], "%li", &loops);
	}

	if (num == tot) {
		for (int i1 = 0; i1 < SNOVA_l * (SNOVA_r - 2); i1++) {
			num_r *= SNOVA_q;
		}
		number_of_keys = 2 * num_r;
	} else {
		// q^{l(r-1)}
		for (int i1 = 0; i1 < SNOVA_l * (SNOVA_r - 1); i1++) {
			num_r *= SNOVA_q;
		}
		start = num * (num_r / tot);

		if (num == (tot - 1)) {
			number_of_keys = num_r - start;
		} else {
			number_of_keys = num_r / tot;
		}
	}

	printf("RectSNOVA (o=%d, r=%d, q=%d, l=%d, m1=%d)  ml2: %d, olr:%d,  Nalpha: %d\n", SNOVA_o, SNOVA_r, SNOVA_q, SNOVA_l,
	       SNOVA_m1, SNOVA_m1 * SNOVA_l2, SNOVA_o * SNOVA_lr, SNOVA_alpha);
	printf("Start: %ld, Tests: %ld / %ld (%.2e), Loops: %ld, Run %ld / %ld\n", start, number_of_keys, num_r, (double)num_r,
	       loops, num, tot);
	fflush(stdout);

	uint64_t cycles0 = 0;
	uint64_t cycles1 = 0;
	uint64_t start_ns = 0;
	uint64_t end_ns = 0;

	for (int i1 = 0; i1 < SNOVA_q; i1++) {
		uint16_t val = 1;
		for (int j1 = 0; j1 < SNOVA_q - 2; j1++) {
			val = (i1 * val) % SNOVA_q;
		}
		invtab[i1] = val;
		mininvtab[i1] = (SNOVA_q - val) % SNOVA_q;
	}

	for (int lind = 0; lind < loops; lind++) {
		uint64_t counts[SNOVA_o * SNOVA_lr + 1] = {0};

		char seed[32];
		if (lind) {
			snprintf((char*)seed, 32, "SNOVA_%ld", (uint64_t)lind);
		} else {
			snprintf((char*)seed, 32, "SNOVA_ABQ");
		}

		uint8_t rand_data[SNOVA_o * SNOVA_lr * EMAT_COLS];
		shake256(rand_data, SNOVA_o * SNOVA_lr * EMAT_COLS, (unsigned char*)seed, strlen(seed));

		uint16_t emat[EMAT_COLS][SNOVA_olr] = {0};

#if 0
		// Random emat
		printf("   *** Random, seed: %s\n", seed);
		for (int i1 = 0; i1 < SNOVA_o * SNOVA_lr; i1++) {
			for (int j1 = 0; j1 < EMAT_COLS; j1++) {
				emat[j1][i1] = rand_data[i1 * EMAT_COLS + j1] % SNOVA_q;
			}
		}
#elif 0
		// Companion emat with random permutation of columns

		uint32_t idx_table[EMAT_COLS] = {0};

		// Permute the Companion
#if 1
		make_idx_table(idx_table, seed);
#else
		for (int mi = 0; mi < EMAT_COLS; ++mi) {
			idx_table[mi] = (229 * mi) % EMAT_COLS;
		}
#endif

#if 0
		printf("uint32_t idx_table[%d] = \n{", EMAT_COLS);
		for (int mi = 0; mi < EMAT_COLS; ++mi) {
			printf("%d, ", idx_table[mi]);
		}
		printf("};\n");
#endif

		uint8_t cmat[SNOVA_o * SNOVA_lr][SNOVA_o * SNOVA_lr] = {0};
		// Companion
		for (int mi = 1; mi < SNOVA_o * SNOVA_lr; ++mi) {
			cmat[mi][mi - 1] = 1;
		}
		cmat[0][SNOVA_o * SNOVA_lr - 1] = 1;

		// Irreducible Polynomials
#if SNOVA_o * SNOVA_lr == 64
		// x^64 + 18*x^7 + 18
		cmat[7][SNOVA_o * SNOVA_lr - 1] = 1;

#elif SNOVA_o * SNOVA_lr == 70
		// x^70 + 18*x^9 + 18*x^4 + 18*x + 18
		cmat[1][SNOVA_o * SNOVA_lr - 1] = 1;
		cmat[4][SNOVA_o * SNOVA_lr - 1] = 1;
		cmat[9][SNOVA_o * SNOVA_lr - 1] = 1;

#elif SNOVA_o * SNOVA_lr == 72
		// x^72 + 18*x^7 + 18*x^5 + 18*x + 18
		cmat[1][SNOVA_o * SNOVA_lr - 1] = 1;
		cmat[5][SNOVA_o * SNOVA_lr - 1] = 1;
		cmat[7][SNOVA_o * SNOVA_lr - 1] = 1;

#else
#error
#endif

		printf("   *** Companion matrix %d, %d, %s\n", SNOVA_o * SNOVA_lr, SNOVA_q, seed);

		// Identity
		emat[idx_table[0]][0] = 1;

		// Derived
		for (int i1 = 1; i1 < EMAT_COLS; ++i1) {
			for (int mi = 0; mi < SNOVA_o * SNOVA_lr; ++mi) {
				uint16_t res = 0;
				for (int mj = 0; mj < SNOVA_o * SNOVA_lr; ++mj) {
					res += cmat[mi][mj] * emat[idx_table[i1 - 1]][mj];
				}
				emat[idx_table[i1]][mi] = res % SNOVA_q;
			}
		}

#if 0
		for (int i1 = 0; i1 < 1 * EMAT_COLS; ++i1) {
			for (int mi = 0; mi < SNOVA_o * SNOVA_lr; ++mi) {
				printf("%3d ", emat[idx_table[i1]][mi]);
			}
			printf("\n");
		}

#endif

#else
		// SNOVA 2.1, ring equation
		printf("   *** Ring E, seed: %s\n", seed);

		uint16_t A[SNOVA_o * SNOVA_alpha * SNOVA_r2];
		uint16_t B[SNOVA_o * SNOVA_alpha * SNOVA_lr];
		uint16_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
		uint16_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
		uint16_t q1[SNOVA_o * SNOVA_alpha * SNOVA_l];
		uint16_t q2[SNOVA_o * SNOVA_alpha * SNOVA_l];

		int index = 0;
		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++)
			for (int i1 = 0; i1 < SNOVA_r2; i1++) {
				A[idx * SNOVA_r2 + i1] = rand_data[index] % SNOVA_q;
				index++;
			}
		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++)
			for (int i1 = 0; i1 < SNOVA_lr; i1++) {
				B[idx * SNOVA_lr + i1] = rand_data[index] % SNOVA_q;
				index++;
			}
		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++)
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				q1[idx * SNOVA_l + i1] = rand_data[index] % SNOVA_q;
				index++;
			}
		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++)
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				q2[idx * SNOVA_l + i1] = rand_data[index] % SNOVA_q;
				index++;
			}

		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
			be_invertible_by_add_aS(&(Am[idx * SNOVA_r2]), &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
			be_invertible_by_add_aS(&(Bm[idx * SNOVA_lr]), &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);

			if (!q1[idx * SNOVA_l + SNOVA_l - 1]) {
				q1[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q1[idx * SNOVA_l] + (q1[idx * SNOVA_l] == 0));
			}
			if (!q2[idx * SNOVA_l + SNOVA_l - 1]) {
				q2[idx * SNOVA_l + SNOVA_l - 1] = SNOVA_q - (q2[idx * SNOVA_l] + (q2[idx * SNOVA_l] == 0));
			}
		}

#if 0
		// Compare to reference implementation
		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
			printf("Am %d\n", idx);
			for (int i1 = 0; i1 < SNOVA_r; i1++) {
				for (int j1 = 0; j1 < SNOVA_r; j1++) {
					printf("%4d", Am[idx * SNOVA_r2 + i1 * SNOVA_r + j1]);
				}
				printf("\n");
			}
			printf("\n");
		}
		printf("\n");

		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
			printf("Bm %d\n", idx);
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				for (int j1 = 0; j1 < SNOVA_r; j1++) {
					printf("%4d", Bm[idx * SNOVA_lr + i1 * SNOVA_r + j1]);
				}
				printf("\n");
			}
			printf("\n");
		}
		printf("\n");

		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
			printf("q1 %d\n", idx);
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				printf("%4d", q1[idx * SNOVA_l + i1]);
			}
			printf("\n");
		}
		printf("\n");

		for (int idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
			printf("q2 %d\n", idx);
			for (int i1 = 0; i1 < SNOVA_l; i1++) {
				printf("%4d", q2[idx * SNOVA_l + i1]);
			}
			printf("\n");
		}
		printf("\n");

#endif

		for (int mi = 0; mi < SNOVA_o; mi++)
			for (int alpha = 0; alpha < SNOVA_alpha; ++alpha) {
				int mj = (alpha + mi) % SNOVA_m1;

				for (int i1 = 0; i1 < SNOVA_r; i1++)
					for (int j1 = 0; j1 < SNOVA_l; j1++)
						for (int a = 0; a < SNOVA_l; a++)
							for (int b = 0; b < SNOVA_l; b++)
								for (int i2 = 0; i2 < SNOVA_r; i2++)
									for (int j2 = 0; j2 < SNOVA_r; j2++) {
										int r1 = (mi * SNOVA_r + i1) * SNOVA_l + j1;
										int r2 = (((mj * SNOVA_r + i2) * SNOVA_l + a) * SNOVA_r + j2) * SNOVA_l + b;

										uint16_t val = q1[(mi * SNOVA_alpha + alpha) * SNOVA_l + a];
										val = (val * q2[(mi * SNOVA_alpha + alpha) * SNOVA_l + b]) % SNOVA_q;
										val = (val * Am[(mi * SNOVA_alpha + alpha) * SNOVA_r2 + i1 * SNOVA_r + i2]) % SNOVA_q;
										val = (val * Bm[(mi * SNOVA_alpha + alpha) * SNOVA_lr + j2 * SNOVA_l + j1]) % SNOVA_q;
										emat[r2][r1] = (emat[r2][r1] + val) % SNOVA_q;
									}
			}

#endif

		// Test coefficients

		for (uint64_t r1 = 0; r1 < number_of_keys; r1++) {
			if (((r1) % 100000000) == 0) {
				if (r1) {
					printf("%ld\n", r1);
				}
				fflush(stdout);
			}

			uint16_t rcoefs[SNOVA_lr] = {0};
			if (num == tot) {
				if (r1 == 0) {
					continue;
				}
				rcoefs[SNOVA_l] = r1 % 2;
				uint64_t r2 = r1 / 2;
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 2; j1 < SNOVA_r; j1++) {
						rcoefs[j1 * SNOVA_l + i1] = r2 % SNOVA_q;
						r2 = r2 / SNOVA_q;
					}
			} else {
				rcoefs[0] = 1;

				uint64_t r2 = start + r1;
				for (int i1 = 0; i1 < SNOVA_l; i1++)
					for (int j1 = 1; j1 < SNOVA_r; j1++) {
						rcoefs[j1 * SNOVA_l + i1] = r2 % SNOVA_q;
						r2 = r2 / SNOVA_q;
					}
			}

			uint16_t etilde[SNOVA_m1l2][SNOVA_olr] = {0};

			start_ns = get_nsec();

			get_etilde(etilde, rcoefs, emat);

			end_ns = get_nsec();
			cycles0 += end_ns - start_ns;
			start_ns = get_nsec();

			int drop = SNOVA_o * SNOVA_lr - get_rank(etilde);

			end_ns = get_nsec();
			cycles1 += end_ns - start_ns;

			counts[drop]++;
		}

		for (int i1 = 0; i1 < 10; i1++) {
			printf("%3d , %8ld\n", i1, counts[i1]);
		}

		uint64_t sum = 0;
		for (int i1 = 10; i1 <= SNOVA_o * SNOVA_lr; i1++) {
			sum += counts[i1];
		}
		printf("10+ , %8ld\n", sum);
	}

	printf("Timing %.3f / %.3f  μsec: %.3f (%.2e)\n", cycles0 / 1e9, cycles1 / 1e9,
	       (cycles0 + cycles1) / 1e3 / number_of_keys / loops, (double)num_r);

	return 0;
}
