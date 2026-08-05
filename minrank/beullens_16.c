// SPDX-License-Identifier: MIT

/**
 * Analyze rank and minrank of E matrix
 *
 * Copyright (c) 2026 SNOVA TEAM
 */

#include <stdalign.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fips202.h"

#define ABQ_ALG1 0
#define ABQ_ALG2 1

#define SNOVA_o 5
#define SNOVA_q 16
#define SNOVA_l 3
#define SNOVA_r 3

#ifndef SNOVA_m1
#define SNOVA_m1 ((SNOVA_o * SNOVA_r) / SNOVA_l)
#endif
#define SNOVA_m2 (SNOVA_o * SNOVA_l * SNOVA_r)

#ifndef SNOVA_alpha
#define SNOVA_alpha (SNOVA_l * SNOVA_r + 2 * SNOVA_r)
#endif

#define SNOVA_l2 (SNOVA_l * SNOVA_l)
#define SNOVA_r2 (SNOVA_r * SNOVA_r)
#define SNOVA_lr (SNOVA_l * SNOVA_r)
#define SNOVA_olr16 ((SNOVA_o * SNOVA_lr + 15) / 16)
#define SNOVA_olr (SNOVA_olr16 * 16)

#define SNOVA_m1l2 (((SNOVA_m1 * SNOVA_l2 + 15) / 16) * 16)

#define SNOVA_ml2 (SNOVA_m1 * SNOVA_l2)

#define EMAT_COLS (SNOVA_m1 * SNOVA_r * SNOVA_l * SNOVA_r * SNOVA_l)

/**
 * Powers of companion matrix to characteristic polynomial of S.
 */
#if SNOVA_l == 2
uint8_t C[2][2][2] = {{{1, 0}, {0, 1}}, {{0, 3}, {1, 14}}};
#elif SNOVA_l == 3
uint8_t C[3][3][3] = {
	{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{0, 0, 14}, {1, 0, 11}, {0, 1, 10}}, {{0, 14, 6}, {0, 11, 12}, {1, 10, 3}}
};
#elif SNOVA_l == 4
uint8_t C[4][4][4] = {{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}},
	{{0, 0, 0, 2}, {1, 0, 0, 11}, {0, 1, 0, 8}, {0, 0, 1, 8}},
	{{0, 0, 2, 3}, {0, 0, 11, 5}, {1, 0, 8, 7}, {0, 1, 8, 4}},
	{{0, 2, 3, 8}, {0, 11, 5, 9}, {0, 8, 7, 3}, {1, 8, 4, 1}}
};
#elif SNOVA_l == 5
uint8_t C[5][5][5] = {{{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}},
	{{0, 0, 0, 0, 12}, {1, 0, 0, 0, 6}, {0, 1, 0, 0, 3}, {0, 0, 1, 0, 15}, {0, 0, 0, 1, 1}},
	{{0, 0, 0, 12, 12}, {0, 0, 0, 6, 10}, {1, 0, 0, 3, 5}, {0, 1, 0, 15, 12}, {0, 0, 1, 1, 14}},
	{{0, 0, 12, 12, 4}, {0, 0, 6, 10, 14}, {0, 0, 3, 5, 11}, {1, 0, 15, 12, 0}, {0, 1, 1, 14, 2}},
	{{0, 12, 12, 4, 11}, {0, 6, 10, 14, 8}, {0, 3, 5, 11, 8}, {0, 15, 12, 0, 6}, {1, 1, 14, 2, 2}}
};
#else
#error "Unsupported parameters"
#endif

uint64_t get_nsec(void) {
	struct timespec time = {0};
	timespec_get(&time, TIME_UTC);
	return (uint64_t)(time.tv_sec * 1e9 + time.tv_nsec);
}

typedef uint8_t gf_t;

static inline uint16_t gf16_expand(const gf_t a) {
	uint16_t val = a | (a << 3) | (a << 6) | (a << 9);
	return val & 0x1111;
}

static inline uint16_t gf16_compress(const uint16_t a) {
	uint16_t val = (a & 0xf) ^ ((a & 0xf0) >> 3) ^ ((a & 0xf00) >> 6) ^ ((a & 0xf000) >> 9);
	return (val ^ ((val & 0xf0) >> 3) ^ (val >> 4)) & 0xf;
}

static inline uint8_t gfni_cleanup(const uint8_t val) {
	return (val ^ ((val & 0xf0) >> 3) ^ (val >> 4)) & 0xf;
}

/**
 * Constant time function. CT is according to valgrind
 */
static inline uint32_t ct_is_not_zero(uint8_t val) {
	// return (val | (val >> 1) | (val >> 2) | (val >> 3)) & 1;
	return val != 0;
}

/**
 * Constant time GF(16) inverse
 *
 * Use that x^q = x and therefore x^(q-2) = x^-1
 */
static inline uint16_t ct_gf_inverse(uint16_t val) {
	uint16_t fact = gf16_compress(val * gf16_expand(val));
	uint16_t res = fact;

	fact = gf16_compress(fact * gf16_expand(fact));
	res = gf16_compress(res * gf16_expand(fact));

	fact = gf16_compress(fact * gf16_expand(fact));
	res = gf16_compress(res * gf16_expand(fact));

	return gf16_expand(res);
}

static gf_t gf_multtab[SNOVA_q * SNOVA_q] = {0};
static gf_t gf_invtab[SNOVA_q] = {0};
static gf_t gf_addtab[SNOVA_q * SNOVA_q] = {0};
static gf_t gf_S[SNOVA_l * SNOVA_l2] = {0};
static uint16_t gf_Sx[SNOVA_l * SNOVA_l2] = {0};

static inline gf_t gf_mult(const gf_t a, const gf_t b) {
	return gf_multtab[a * SNOVA_q + b];
}

static inline gf_t gf_inv(const gf_t a) {
	return gf_invtab[a];
}

static inline gf_t gf_add(const gf_t a, const gf_t b) {
	return gf_addtab[a * SNOVA_q + b];
}

static inline void gf_set_add(gf_t* a, const gf_t b) {
	*a = gf_addtab[*a * SNOVA_q + b];
}

static inline gf_t gf_sub(const gf_t a, const gf_t b) {
#if SNOVA_q != 16
	return gf_addtab[a * SNOVA_q + (SNOVA_q - b) % SNOVA_q];
#else
	return gf_addtab[a * SNOVA_q + b];
#endif
}

static inline void init_gf_tables(void) {
	// GF(16)
	uint8_t F_star[15] = {1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9};  // Z2[x]/(x^4+x+1)
	for (int i1 = 0; i1 < 16; i1++) {
		gf_multtab[i1] = 0;
		gf_multtab[i1 * SNOVA_q] = 0;
	}
	for (int i1 = 0; i1 < SNOVA_q - 1; i1++)
		for (int j1 = 0; j1 < SNOVA_q - 1; j1++) {
			gf_multtab[F_star[i1] * SNOVA_q + F_star[j1]] = F_star[(i1 + j1) % (SNOVA_q - 1)];
		}

	for (int i1 = 0; i1 < SNOVA_q; i1++)
		for (int j1 = 0; j1 < SNOVA_q; j1++) {
			gf_addtab[i1 * SNOVA_q + j1] = (i1 ^ j1);
		}
	// Use that x^q = x and therefore x^(q-2) = x^-1
	for (int i1 = 0; i1 < SNOVA_q; i1++) {
		gf_t val = i1;
		for (int j1 = 3; j1 < SNOVA_q; j1++) {
			val = gf_mult(val, i1);
		}
		gf_invtab[i1] = val;
	}
}

static inline void gf_mat_mul(gf_t* a, const gf_t* b, const gf_t* c) {
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_t sum = 0;
			for (int k1 = 0; k1 < SNOVA_l; k1++) {
				gf_set_add(&sum, gf_mult(b[i1 * SNOVA_l + k1], c[k1 * SNOVA_l + j1]));
			}
			a[i1 * SNOVA_l + j1] = sum;
		}
}

// Set the irreducible S matrix
static inline void set_S(gf_t* gf_S1) {
	for (int i1 = 0; i1 < SNOVA_l; i1++)
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_S1[i1 * SNOVA_l + j1] = 8 - (i1 + j1);
		}
#if SNOVA_l == 5
	gf_S1[SNOVA_l2 - 1] = 9;
#endif
}

static inline void gen_S_array(void) {
	memset(gf_S, 0, sizeof(gf_S));

	for (int i1 = 0; i1 < SNOVA_l; i1++) {
		gf_S[i1 * SNOVA_l + i1] = 1;
	}

#if SNOVA_l > 1
	set_S(&gf_S[1 * SNOVA_l2]);

	for (int i1 = 2; i1 < SNOVA_l; i1++) {
		gf_mat_mul(&gf_S[i1 * SNOVA_l2], &gf_S[1 * SNOVA_l2], &gf_S[(i1 - 1) * SNOVA_l2]);
	}
#endif

	for (int i1 = 0; i1 < SNOVA_l * SNOVA_l2; i1++) {
		gf_Sx[i1] = gf16_expand(gf_S[i1]);
	}
}

/**
 * Utilities
 */

#if SNOVA_l == 4

static inline gf_t gf_mat_det(gf_t* a) {
#define DET_SUB(a, b) (a ^ b)
#define DET_MULT(a, b) gf_multtab[a * SNOVA_q + b]
	gf_t det = 0;
	gf_t det_l;
	gf_t det_r;
#define DET_L(x, y) det_l = DET_SUB(DET_MULT(a[x], a[4 + y]), DET_MULT(a[y], a[4 + x]))
#define DET_R(x, y) det_r = DET_SUB(DET_MULT(a[8 + x], a[12 + y]), DET_MULT(a[8 + y], a[12 + x]))
#define DET22(x1, y1, x2, y2) \
    DET_L(x1, y1);            \
    DET_R(x2, y2);            \
    det ^= DET_MULT(det_l, det_r)
	DET22(0, 1, 2, 3);
	DET22(0, 2, 3, 1);
	DET22(0, 3, 1, 2);
	DET22(1, 2, 0, 3);
	DET22(1, 3, 2, 0);
	DET22(2, 3, 0, 1);

	return det;
}

#else
static inline gf_t gf_mat_det(gf_t* a) {
	gf_t det = 0;
#if SNOVA_l == 1
	det = a[0];
#elif SNOVA_l == 2
	det = gf_sub(gf_mult(a[0], a[3]), gf_mult(a[1], a[2]));
#elif SNOVA_l == 3
	det = gf_mult(a[0], gf_sub(gf_mult(a[4], a[8]), gf_mult(a[5], a[7])));
	gf_set_add(&det, gf_mult(a[1], gf_sub(gf_mult(a[5], a[6]), gf_mult(a[3], a[8]))));
	gf_set_add(&det, gf_mult(a[2], gf_sub(gf_mult(a[3], a[7]), gf_mult(a[4], a[6]))));
#elif SNOVA_l == 5
	gf_t det_l;
	gf_t det_r;
#define DET_L(x, y) det_l = gf_sub(gf_mult(a[x], a[5 + y]), gf_mult(a[y], a[5 + x]))
#define DET_R2(x, y, z) gf_mult(gf_sub(gf_mult(a[10 + x], a[15 + y]), gf_mult(a[10 + y], a[15 + x])), a[20 + z])
#define DET_R3(x, y, z) det_r = gf_add(DET_R2(x, y, z), gf_add(DET_R2(y, z, x), DET_R2(z, x, y)))
#define DET23(x1, y1, x2, y2, z2) \
    DET_L(x1, y1);                \
    DET_R3(x2, y2, z2);           \
    gf_set_add(&det, gf_mult(det_l, det_r))
	DET23(0, 1, 2, 3, 4);
	DET23(0, 2, 3, 1, 4);
	DET23(0, 3, 1, 2, 4);
	DET23(0, 4, 1, 3, 2);
	DET23(1, 2, 0, 3, 4);
	DET23(1, 3, 2, 0, 4);
	DET23(1, 4, 2, 3, 0);
	DET23(2, 3, 0, 1, 4);
	DET23(2, 4, 0, 3, 1);
	DET23(3, 4, 2, 0, 1);
#undef DET_R2
#undef DET_R3
#undef DET23
#undef DET_L
#else
#error "Unsupported rank"
#endif
	return det;
}
#endif

void convert_bytes_to_GF16s(const uint8_t* byte_array, gf_t* gf16_array, int num_of_GF16s) {
	int i;
	int pairs = num_of_GF16s >> 1;

	// Convert each byte into two GF16 values
	for (i = 0; i < pairs; ++i) {
		gf16_array[i * 2] = byte_array[i] & 0x0F;
		gf16_array[i * 2 + 1] = (byte_array[i] >> 4) & 0x0F;
	}

	// Handle the last GF16 value if num_of_GF16s is odd
	if (num_of_GF16s % 2 == 1) {
		gf16_array[num_of_GF16s - 1] = byte_array[pairs] & 0x0F;
	}
}

uint8_t fixed_abq[SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l)] = {0};

static int gen_fixed_ABQ(const char* abq_seed) {
	uint8_t rng_out[SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l)] = {0};

	shake256(rng_out, SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l), (uint8_t*)abq_seed, strlen(abq_seed));
	convert_bytes_to_GF16s(rng_out, fixed_abq, SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr + 2 * SNOVA_l));

	int fails = 0;
#if !ABQ_ALG1
	// Check if q1 and q2 are always non zero
	gf_t* aptr = fixed_abq;
	gf_t* q1 = aptr + SNOVA_o * SNOVA_alpha * (SNOVA_r2 + SNOVA_lr);
	gf_t* q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	// Check if q1 and q2 are always non zero
	for (int i1 = 0 ; i1 < SNOVA_o * SNOVA_alpha; i1++) {
		uint32_t sum1 = 0;
		uint32_t sum2 = 0;
		for (int j1 = 0 ; j1 < SNOVA_l; j1 ++) {
			sum1 |= q1[i1 * SNOVA_l + j1];
			sum2 |= q2[i1 * SNOVA_l + j1];
		}
		if ((sum1 == 0) || (sum2 == 0)) {
			fails++;
		}
	}
	if (fails) {
		static int first = 1;
		if (first) {
			printf("Warning: Some q1,q2 are zero for (%s): %d\n", abq_seed, fails);
			first = 0;
		}
	}
#endif
	return fails;
}

static inline void snova_init(void) {
	init_gf_tables();
	gen_S_array();
}

/**
 * Ensure that a matrix is invertible by adding multiples of S
 */
static inline void be_invertible_by_add_aS(gf_t* mat, const gf_t* orig, const int l1, const int l2) {
	memcpy(mat, orig, l1 * l2);
#if ABQ_ALG2 && SNOVA_l > 1
	if ((l1 == SNOVA_l) && (l2 == SNOVA_l))
		if (gf_mat_det(mat) == 0) {
			for (gf_t f1 = 1; f1 < SNOVA_q; f1++) {
				for (int i1 = 0; i1 < SNOVA_l2; i1++) {
					gf_set_add(&mat[i1], gf_mult(f1, gf_S[SNOVA_l2 + i1]));
				}
				if (gf_mat_det(mat) != 0) {
					break;
				}
			}
		}
#endif
}

/**
 * Improve q and calculate Q matrix
 */
static inline void gen_a_FqS(gf_t* Qm, gf_t* q) {
#if ABQ_ALG1
	if (!q[SNOVA_l - 1]) {
		q[SNOVA_l - 1] = SNOVA_q - (q[0] + (q[0] == 0));
	}
#endif

	for (int i1 = 0; i1 < SNOVA_l2; i1++) {
		gf_t sum = 0;
		for (int j1 = 0; j1 < SNOVA_l; j1++) {
			gf_set_add(&sum, gf_mult(q[j1], gf_S[j1 * SNOVA_l2 + i1]));
		}
		Qm[i1] = sum;
	}
}

/**
 * Use last part of the P matrix to establish ABQ
 */
static int gen_ABQ(char* seed, gf_t* A, gf_t* Am, gf_t* Bm, gf_t* Q1m, gf_t* Q2m) {
	int res = gen_fixed_ABQ(seed);

	gf_t* B = A + SNOVA_o * SNOVA_alpha * SNOVA_r2;
	gf_t* q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
	gf_t* q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

	for (size_t idx = 0; idx < SNOVA_o * SNOVA_alpha; idx++) {
		be_invertible_by_add_aS(&Am[idx * SNOVA_r2], &A[idx * SNOVA_r2], SNOVA_r, SNOVA_r);
		be_invertible_by_add_aS(&Bm[idx * SNOVA_lr], &B[idx * SNOVA_lr], SNOVA_r, SNOVA_l);
		gen_a_FqS(&Q1m[idx * SNOVA_l2], &q1[idx * SNOVA_l]);
		gen_a_FqS(&Q2m[idx * SNOVA_l2], &q2[idx * SNOVA_l]);
	}
	return res;
}

// Generate ABQ
gf_t Am[SNOVA_o * SNOVA_alpha * SNOVA_r2];
gf_t Bm[SNOVA_o * SNOVA_alpha * SNOVA_lr];
gf_t Q1[SNOVA_o * SNOVA_alpha * SNOVA_l2];
gf_t Q2[SNOVA_o * SNOVA_alpha * SNOVA_l2];

gf_t *q1;
gf_t *q2;
uint16_t z1bx[SNOVA_o * SNOVA_alpha * SNOVA_l2] = {0};
uint16_t z2bx[SNOVA_o * SNOVA_alpha * SNOVA_l2] = {0};

void get_etilde(uint16_t etilde[SNOVA_m1l2][SNOVA_olr], uint16_t rcoefs[SNOVA_lr]) {
	memset(etilde, 0, SNOVA_m1l2 * SNOVA_olr * sizeof(uint16_t));

	alignas(32) uint16_t etildex[SNOVA_ml2 * SNOVA_olr] = {0};

	for (int mi = 0; mi < SNOVA_o; mi++)
		for (int alpha = 0; alpha < SNOVA_alpha; alpha++) {
			int mj = (alpha + mi) % SNOVA_m1;
			int mia = mi * SNOVA_alpha + alpha;

			alignas(32) uint16_t z1a[SNOVA_lr] = {0};
			for (int i1 = 0; i1 < SNOVA_r; i1++)
				for (int a1 = 0; a1 < SNOVA_l; a1++)
					for (int i2 = 0; i2 < SNOVA_r; i2++) {
						z1a[i1 * SNOVA_l + a1] ^= gf_mult(Am[mia * SNOVA_r2 + i1 * SNOVA_r + i2], rcoefs[i2 * SNOVA_l + a1]);
					}

			alignas(32) uint16_t z1[SNOVA_lr] = {0};
			for (int i1 = 0; i1 < SNOVA_r; i1++)
				for (int a2 = 0; a2 < SNOVA_l; a2++)
					for (int a1 = 0; a1 < SNOVA_l; a1++) {
						z1[i1 * SNOVA_l + a2] ^= gf_mult(z1a[i1 * SNOVA_l + a1], z1bx[mia * SNOVA_l2 + a1 * SNOVA_l + a2]);
					}

			alignas(32) uint16_t z2a[SNOVA_l2] = {0};
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int a1 = 0; a1 < SNOVA_l; a1++)
					for (int i2 = 0; i2 < SNOVA_r; i2++) {
						z2a[i1 * SNOVA_l + a1] ^= gf_mult(Bm[mia * SNOVA_lr + i2 * SNOVA_l + i1], rcoefs[i2 * SNOVA_l + a1]);
					}

			alignas(32) uint16_t z2[SNOVA_l2] = {0};
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int a2 = 0; a2 < SNOVA_l; a2++)
					for (int a1 = 0; a1 < SNOVA_l; a1++) {
						z2[i1 * SNOVA_l + a2] ^= gf_mult(z2a[i1 * SNOVA_l + a1], z2bx[mia * SNOVA_l2 + a1 * SNOVA_l + a2]);
					}

			alignas(32) uint16_t z1x[SNOVA_lr] = {0};
			for (int i1 = 0; i1 < SNOVA_r; i1++)
				for (int a2 = 0; a2 < SNOVA_l; a2++) {
					z1x[i1 * SNOVA_l + a2] = gf16_expand(z1[i1 * SNOVA_l + a2]);
				}

			for (int i1 = 0; i1 < SNOVA_r; i1++)
				for (int a2 = 0; a2 < SNOVA_l; a2++)
					for (int j1 = 0; j1 < SNOVA_l2; j1++)
						etildex[((mi * SNOVA_r + i1) * SNOVA_l + a2) * SNOVA_ml2 + mj * SNOVA_l2 + j1] ^=
						    z1x[i1 * SNOVA_l + a2] * z2[j1];
		}

	for (int mi = 0; mi < SNOVA_o * SNOVA_r; mi++)
		for (int mj = 0; mj < SNOVA_m1; mj++)
			for (int a2 = 0; a2 < SNOVA_l; a2++)
				for (int j1 = 0; j1 < SNOVA_l; j1++)
					for (int b2 = 0; b2 < SNOVA_l; b2++)
						etilde[(mj * SNOVA_l + a2) * SNOVA_l + b2][mi * SNOVA_l + j1] =
						    gf16_compress(etildex[(mi * SNOVA_l + a2) * SNOVA_ml2 + mj * SNOVA_l2 + j1 * SNOVA_l + b2]);
}

int get_rank(uint16_t matrix[SNOVA_m1l2][SNOVA_olr]) {
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

		uint16_t t_GF16 = ct_gf_inverse(matrix[rank][j1]);

		for (int i2 = rank + 1; i2 < SNOVA_m1 * SNOVA_l2; i2++) {
			uint16_t gji = gf16_expand(gf16_compress(t_GF16 * matrix[i2][j1]));
			for (int j2 = 0; j2 < SNOVA_olr; j2++) {
				matrix[i2][j2] ^= gji * matrix[rank][j2];
			}
		}

		for (int i2 = rank + 1; i2 < SNOVA_m1 * SNOVA_l2; i2++) {
			for (int j2 = 0; j2 < SNOVA_olr; j2++) {
				matrix[i2][j2] = gf16_compress(matrix[i2][j2]);
			}
		}

		rank++;
	}

	return rank;
}

int get_rank_tr(uint16_t matrix[SNOVA_olr][SNOVA_m1l2]) {
	int rank = 0;
	for (int j1 = 0; j1 < SNOVA_m1 * SNOVA_l2; j1++) {
		int i1;

		for (i1 = rank; i1 < SNOVA_o * SNOVA_l * SNOVA_r - 1; i1++) {
			if (matrix[i1][j1]) {
				break;
			}
		}

		if (matrix[i1][j1] == 0) {
			continue;
		}

		if (i1 > rank) {
			for (int j2 = 0; j2 < SNOVA_m1l2; j2++) {
				uint16_t temp = matrix[rank][j2];
				matrix[rank][j2] = matrix[i1][j2];
				matrix[i1][j2] = temp;
			}
		}

		uint16_t t_GF16 = ct_gf_inverse(matrix[rank][j1]);

		for (int i2 = rank + 1; i2 < SNOVA_o * SNOVA_l * SNOVA_r; i2++) {
			uint16_t gji = gf16_expand(gf16_compress(t_GF16 * matrix[i2][j1]));
			for (int j2 = 0; j2 < SNOVA_m1l2; j2++) {
				matrix[i2][j2] ^= gji * matrix[rank][j2];
			}
		}

		for (int i2 = rank + 1; i2 < SNOVA_o * SNOVA_l * SNOVA_r; i2++) {
			for (int j2 = 0; j2 < SNOVA_m1l2; j2++) {
				matrix[i2][j2] = gf16_compress(matrix[i2][j2]);
			}
		}

		rank++;
	}

	return rank;
}

int main(int argc, char** argv) {
	if (argc < 3 || argc > 4) {
		printf("Supply index, partitions, and (optional) number of loops\n");
		printf("\tUsing an index of -1 will analyze the scalar subset\n");
		printf("\tUsing an index of -2 will analyze the remaining subset with a_00 = 0\n");
		return 0;
	}

	double num_r = 1;
	uint64_t start = 0;
	uint64_t number_of_keys;

	int64_t num = 0;
	int64_t tot = 1;
	uint64_t loops = 1;

	sscanf(argv[1], "%li", &num);
	sscanf(argv[2], "%li", &tot);
	if (num >= tot) {
		printf("Error: index \u2265 number of partitions\n");
		return 0;
	}
	if (argc == 4) {
		sscanf(argv[3], "%lu", &loops);
	}

	if (num < 0) {
		for (int i1 = 0; i1 < SNOVA_l * (SNOVA_r - 1); i1++) {
			num_r *= SNOVA_q;
		}
		number_of_keys = 1;
		for (int i1 = 0; i1 < (SNOVA_r - 1); i1++) {
			number_of_keys *= SNOVA_q;
		}
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

	printf("RectSNOVA (o=%d, q=%d, l=%d, r=%d, m1=%d)  ml2: %d, olr:%d,  Nalpha: %d\n", SNOVA_o, SNOVA_q, SNOVA_l, SNOVA_r,
	       SNOVA_m1, SNOVA_m1 * SNOVA_l2, SNOVA_o * SNOVA_lr, SNOVA_alpha);
	printf("Start: %ld, Tests: %ld / %.0f (%.2e), Loops: %ld, Run %ld / %ld\n", start, number_of_keys, num_r, num_r, loops, num,
	       tot);
	printf("Alg 1: %d, and alg 2: %d\n", ABQ_ALG1, ABQ_ALG2);
	fflush(stdout);

	snova_init();

	uint64_t cycles0 = 0;
	uint64_t cycles1 = 0;
	uint64_t cycles2 = 0;

	uint64_t ranks[SNOVA_m2 + 1] = {0};
	int minminrank = SNOVA_m2 + 1;
	int maxrank = SNOVA_m2;
	if (SNOVA_m1 * SNOVA_l2 < SNOVA_m2) {
		maxrank = SNOVA_m1 * SNOVA_l2;
	}

	char seed[32] = {0};
	uint64_t lindex = 0;
	for (uint64_t index = 0; index < loops; index++) {
		if (index) {
			snprintf((char*)seed, 32, "SNOVA_ABQ_%ld", (uint64_t)lindex);
		} else {
			snprintf((char*)seed, 32, "SNOVA_ABQ");
		}
		if (index) {
			// printf("Using %s\n", seed);
		}
		lindex++;

		if (gen_ABQ(seed, fixed_abq, Am, Bm, Q1, Q2)) {
			// printf("Skipping %s\n", seed);
			index--;
			continue;
		};

		gf_t* B = fixed_abq + SNOVA_o * SNOVA_alpha * SNOVA_r2;
		q1 = B + SNOVA_o * SNOVA_alpha * SNOVA_lr;
		q2 = q1 + SNOVA_o * SNOVA_alpha * SNOVA_l;

		memset(z1bx, 0, sizeof(z1bx));
		for (int mia = 0; mia < SNOVA_o * SNOVA_alpha; mia++)
			for (int a1 = 0; a1 < SNOVA_l; a1++)
				for (int a2 = 0; a2 < SNOVA_l; a2++)
					for (int a = 0; a < SNOVA_l; a++) {
						z1bx[mia * SNOVA_l2 + a1 * SNOVA_l + a2] ^= gf_mult(C[a1][a2][a], q1[mia * SNOVA_l + a]);
					}

		memset(z2bx, 0, sizeof(z1bx));
		for (int mia = 0; mia < SNOVA_o * SNOVA_alpha; mia++)
			for (int a1 = 0; a1 < SNOVA_l; a1++)
				for (int a2 = 0; a2 < SNOVA_l; a2++)
					for (int a = 0; a < SNOVA_l; a++) {
						z2bx[mia * SNOVA_l2 + a1 * SNOVA_l + a2] ^= gf_mult(C[a1][a2][a], q2[mia * SNOVA_l + a]);
					}

		for (uint64_t r1 = 0; r1 < number_of_keys; r1++) {
			if (((r1) % 100000000) == 0) {
				if (r1) {
					printf("%ld\n", r1);
				}
				fflush(stdout);
			}

			uint16_t rcoefs[SNOVA_lr] = {0};

			if (num == -2 && r1 == 0) {
				continue;
			}
			if (num != -2) {
				rcoefs[0] = 1;
			}

			uint64_t r2 = start + r1;
			for (int i1 = 0; i1 < SNOVA_l; i1++)
				for (int j1 = 1; j1 < SNOVA_r; j1++) {
					rcoefs[j1 * SNOVA_l + i1] = r2 & 0xf;
					r2 = r2 >> 4;
				}

			uint16_t matrix[SNOVA_m1l2][SNOVA_olr] = {0};

			cycles0 = get_nsec();
			get_etilde(matrix, rcoefs);
			cycles1 += get_nsec() - cycles0;

#if SNOVA_m1 * SNOVA_l2 < SNOVA_o * SNOVA_l * SNOVA_r
			uint16_t matrix_tr[SNOVA_olr][SNOVA_m1l2] = {0};
			for (int i1 = 0; i1 < SNOVA_m1 * SNOVA_l2; i1++)
				for (int j1 = 0; j1 < SNOVA_o * SNOVA_l * SNOVA_r; j1++) {
					matrix_tr[j1][i1] = matrix[i1][j1];
				}

			cycles0 = get_nsec();
			int rank = get_rank_tr(matrix_tr);
			cycles2 += get_nsec() - cycles0;
#else

			cycles0 = get_nsec();
			int rank = get_rank(matrix);
			cycles2 += get_nsec() - cycles0;
#endif

			if (rank > maxrank) {
				printf("Rank Error %d (%ld) '%s' max:%d\n", rank, r1, seed, maxrank);
				break;
			}

			ranks[rank] += 1;
			if (rank < minminrank) {
				minminrank = rank;
			}

			if (rank <= maxrank - 2 * SNOVA_l) {
				printf("Rank %d (%ld) '%s' %ld\n", rank, r1, seed, index);
				fflush(stdout);
				// break;
			}
		}
	}

	for (int idx = maxrank; idx >= minminrank; idx--) {
		printf("%d, %8ld, %.03f\n", maxrank - idx, ranks[idx], (double)ranks[idx] / (double)ranks[maxrank]);
	}

	printf("Timing %.3f / %.3f  μsec: %.3f\n", cycles1 / 1e9, cycles2 / 1e9,
	       (cycles1 + cycles2) / 1e3 / number_of_keys / loops);
	printf("  Skipped %ld seeds\n", lindex - loops);

	return 0;
}
