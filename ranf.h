
#ifndef RANF_H
#define RANF_H

#include <stddef.h>
#include <inttypes.h>

#include "meowhash.h"

typedef int32_t int_ranf_t;
#define MEOW_EXTRACT MeowU32From
typedef int_ranf_t * ranf_state;

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

#define CACHE_LINE_COUNT (CACHE_LINE_SIZE/sizeof(int_ranf_t))

#define CACHE_INDEX(thread) ((thread)*CACHE_LINE_COUNT)

/*
**********************************************************************
     ranf_state create_state_mt(
       const size_t thread_count, const long seed1, const long seed2
     );
     Allocate and initialize random number generator state(s)
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
const ranf_state create_salt_mt(const size_t thread_count, const int_ranf_t seed1, const int_ranf_t seed2);
ranf_state create_seed_mt(const ranf_state salt, const size_t thread_count);

#define hash_mt(STATE, THREAD, SALT, ARGC, ...) {\
  const size_t curntg = CACHE_INDEX(THREAD); \
  int args[ARGC] = { __VA_ARGS__ }; \
  const meow_u128 Hash = MeowHash(SALT + curntg, ARGC * sizeof(int), args); \
  STATE[curntg] = MEOW_EXTRACT(Hash, 0); \
  STATE[curntg + 1] = MEOW_EXTRACT(Hash, 1); \
}

/* RANDLIB functions */

/* draw a random U[0,1] deviate */
/* state: the random number generator state */
/* thread: which thread is accessing the state */
double ranf_mt(ranf_state state, const size_t thread);
/* draw a random U[0,1] deviate from thread 0 state */
double ranf(ranf_state state);

size_t ignbin_mt(ranf_state state, const size_t n, const double pp, const size_t thread);
size_t ignbin(ranf_state state, const size_t n, const double pp);

size_t ignpoi_mt(ranf_state state, const double lambda, const size_t thread);
size_t ignpoi(ranf_state state, const double lambda);

double sexpo_mt(ranf_state state, const size_t tn);
double sexpo(ranf_state state);

double snorm(ranf_state state);
double snorm_mt(ranf_state state, const size_t tn);

double fsign(const double num, const double sign);

double gengam(ranf_state state, const double a, const double r);
double gengam_mt(ranf_state state, const double a, const double r, const size_t tn);

//added some new beta, gamma generating functions: ggilani 27/11/14
double gen_norm(ranf_state state, const double mu, const double sigma);
double gen_norm_mt(ranf_state state, const double mu, const double sigma, const size_t tn);

double gen_gamma(ranf_state state, const double a, const double r);
double gen_gamma_mt(ranf_state state, const double a, const double r, const size_t tn);

double gen_beta(ranf_state state, const double alpha, const double beta);
double gen_beta_mt(ranf_state state, const double alpha, const double beta, const size_t tn);

//added some new lognormal sampling functions: ggilani 09/02/17
double gen_lognormal(double,double);
double gen_lognormal_mt(double,double,int);
double sgamma(double);
double sgamma_mt(double,int);

double gammln(const double xx);

// void SampleWithoutReplacement(int,int,int);
double WeibullCDF(const int x, const double shape, const double scale);

void FYShuffle_mt(ranf_state state, int * array, const size_t n, const size_t tn);

#endif // RANF_H