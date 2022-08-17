/**
 * RNG implementation.
 * Derived from https://svn.r-project.org/R/
 */

#include "math_ext.h"

#include <stdint.h>
#include <math.h>

static unsigned int i_seed[628];                        /* allow for optimizing compilers to read over bound */

#define I1                      (i_seed[0])
#define I2                      (i_seed[1])

#define N                       624
#define M                       397
#define ALPHA                   69069                   /* https://www.gnu.org/software/gsl/doc/html/rng.html */
#define INIT_SCRAMBLING         50
#define UPPER_MASK              0x80000000              /* most significant w-r bits */
#define LOWER_MASK              0x7fffffff              /* least significant r bits */
#define i2_32m1                 2.328306437080797e-10   /* = 1/(2^32 - 1) */
#define MATRIX_A                0x9908b0df              /* constant vector a */
#define DEFAULT_SGEN_SEED       4357

#define TEMPERING_MASK_B        0x9d2c5680
#define TEMPERING_MASK_C        0xefc60000
#define TEMPERING_SHIFT_U(y)    ((y) >> 11)
#define TEMPERING_SHIFT_S(y)    ((y) <<  7)
#define TEMPERING_SHIFT_T(y)    ((y) << 15)
#define TEMPERING_SHIFT_L(y)    ((y) >> 18)

static unsigned int mti = N + 1;                        /* mti==N+1 means mt[N] is not initialized */
static unsigned int *mt = i_seed + 1;                   /* the array for the state vector  */

/**
 * Initialize seed environment using "Mersenne-Twister" approach.
 * @param seed Seed.
 */
static void seed_init(unsigned int seed) {
    /* Initial scrambling */
    for (int j = 0; j < INIT_SCRAMBLING; j++) {
        seed = (ALPHA * seed + 1);
    }

    /* Mersenne-Twister */
    for (int j = 0; j < 1 + N; j++) {
        seed = (ALPHA * seed + 1);
        i_seed[j] = seed;
    }

    I1         = N;
    if (I1 <= 0) { I1 = N; }

    /* check for all zeroes */
    for (int j = 1; j <= N; j++) {
        if (i_seed[j] != 0) { break; }
    }
}

static void sgen_rand(unsigned int seed) {
    for (int i = 0; i < N; i++) {
        mt[i] = seed & 0xffff0000;
        seed = ALPHA * seed + 1;
        mt[i] |= (seed & 0xffff0000) >> 16;
        seed = ALPHA * seed + 1;
    }
    mti        = N;
}

/**
 * @return reals: [0,1)-interval
 */
static double gen_rand() {
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    static unsigned int mag01[2] = {0x0, MATRIX_A};
    unsigned int        y;
    mti = i_seed[0];

    /* generate N words at one time */
    if (mti >= N) {
        int kk;

        /* if sgenrand() has not been called, a default initial seed is used */
        if (mti == N + 1) { sgen_rand(DEFAULT_SGEN_SEED); }

        for (kk = 0; kk < N - M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (; kk < N - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }

    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    i_seed[0] = mti;

    return ((double) y * 2.3283064365386963e-10);
}

double unif_rand() {
    double x = gen_rand();
    if (x <= 0.0) { return 0.5 * i2_32m1; }
    if ((1.0 - x) <= 0.0) { return 1.0 - 0.5 * i2_32m1; }
    return x;
}

void set_seed(const unsigned int seed) {
    seed_init(seed);
}

static double rbits(int bits) {
    uint64_t v = 0;

    for (int n = 0; n <= bits; n += 16) {
        int v1 = (int) floor(unif_rand() * 65536);
        v = 65536 * v + v1;
    }

    const int_least64_t one64 = 1L;
    return (double) (v & ((one64 << bits) - 1));
}

int unif_index(int dn) {
    if (dn <= 0) { return 0; }
    int    bits = (int) ceil(log2(dn));
    double dv;
    do { dv = rbits(bits); } while (dn <= dv);
    return (int) dv;
}

int_arr_2_t unif_index2(int dn) {
    int_arr_2_t ry;
    ry.arr[0] = -1;
    ry.arr[1] = -1;
    if (dn < 2) { return ry; }

    int x[dn], n = dn;

    for (int i = 0; i < dn; i++) { x[i] = i; }
    for (int i = 0; i < 2; i++) {
        int j     = unif_index(n);
        ry.arr[i] = x[j];
        x[j]      = x[--n];
    }
    return ry;
}
