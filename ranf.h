
#ifndef RANF_H
#define RANF_H

typedef long * ranf_state;

long * create_state_mt(int thread_count, long seed1, long seed2);

/* draw a random U[0,1] deviate */
/* state: the random number generator state */
/* thread: which thread is accessing the state */
double ranf_mt(long * state, int thread);
/* draw a random U[0,1] deviate from thread 0 state */
double ranf(long * state);

void hash_mt(long * state, int thread, long * salt, int argc, ...);

/* RANDLIB functions */
long ignbin_mt(long, double, int);
long ignbin(long, double);

long ignpoi(double );

long ignpoi_mt(double,int);

void setall(long ,long);
double sexpo_mt(int);
double sexpo(void);
long mltmod(long ,long ,long );
double snorm(void);
double snorm_mt(int);
double fsign(double, double );
double gengam(double,double);
double gengam_mt(double,double,int);
//added some new beta, gamma generating functions: ggilani 27/11/14
double gen_norm(double,double);
double gen_norm_mt(double,double,int);
double gen_gamma(double,double);
double gen_gamma_mt(double,double,int);
double gen_beta(double, double);
double gen_beta_mt(double,double,int);
//added some new lognormal sampling functions: ggilani 09/02/17
double gen_lognormal(double,double);
double gen_lognormal_mt(double,double,int);
double sgamma(double);
double sgamma_mt(double,int);
void SampleWithoutReplacement(int,int,int);
double gammln(double );
double WeibullCDF(int,double,double);

#endif // RANF_H