
// CACHE_LINE_SIZE should be set to the number of longs that fit in a cache line
// on 64 bit systems (8 byte longs) this is typically 64 bytes, so CACHE_LINE_SIZE=8

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ranf.h"

/* ranf defines */
#define Xm1 2147483563
#define Xm2 2147483399
#define Xa1 40014
#define Xa2 40692
#define Xa1vw 2082007225
#define Xa2vw 84306273

// would be used as initial seed >> hash on the elements of ... => new seed set, then rand_mt
// assuming ...integer inputs, make them into longs, then push them thru as hashs, with initial seed as salt
// void hash_mt(
// 	ranf_state state, const int thread, const ranf_state salt, const size_t argc, ...
// ) {
// 	va_list args;
// 	va_start(args, argc);

// 	const size_t curntg = CACHE_INDEX(thread);
// 	meow_u128 Hash = MeowHash(salt[curntg], Size, Buffer);
//     state[curntg] = MeowU64From(Hash, 0);
// 	state[curntg + 1] = MeowU64From(Hash, 1);



// 	va_end(args);
// }

double ranf_mt(ranf_state state, const size_t thread) {
    const size_t curntg = CACHE_INDEX(thread);

	int_ranf_t s1 = state[curntg];
	int_ranf_t s2 = state[curntg + 1];
	int_ranf_t k = s1 / 53668L;

	s1 = Xa1 * (s1 - k * 53668L) - k * 12211;
	if (s1 < 0) s1 += Xm1;
	k = s2 / 52774L;
	s2 = Xa2 * (s2 - k * 52774L) - k * 3791;
	if (s2 < 0) s2 += Xm2;

	state[curntg] = s1;
	state[curntg + 1] = s2;
	
    long z = s1 - s2;
	if (z < 1) z += (Xm1-1);
	return ((double) z)/Xm1;
}

double ranf(ranf_state state) { return ranf_mt(state, 0); }

const ranf_state create_salt_mt(const size_t thread_count, const int_ranf_t seed1, const int_ranf_t seed2) {

    ranf_state state = malloc(sizeof(int_ranf_t) * thread_count * CACHE_LINE_COUNT);
    state[0] = seed1;
    state[1] = seed2;
    for(size_t g = 1; g < thread_count; g++) {
        state[CACHE_INDEX(g)] = Xa1vw * state[CACHE_INDEX(g-1)] % Xm1;
        state[CACHE_INDEX(g) + 1] = Xa2vw * state[CACHE_INDEX(g-1) + 1] % Xm2;
    }
    return state;
}

ranf_state create_seed_mt(const ranf_state salt, const size_t thread_count) {

    ranf_state state = malloc(sizeof(uint64_t) * thread_count * CACHE_LINE_COUNT);
	memcpy(state, salt, sizeof(uint64_t) * thread_count * CACHE_LINE_COUNT);
    return state;
}

size_t ignbin(ranf_state state, const size_t n, const double pp) {
    return ignbin_mt(state, n, pp, 0);
};

#define ERR_CRITICAL(x) {fprintf(stderr,x);exit(1);}
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define minF(a,b) ((a) <= (b) ? (a) : (b))
#define maxF(a,b) ((a) >= (b) ? (a) : (b))

size_t ignbin_mt(ranf_state state, const size_t n, const double pp, const size_t tn)
/*
**********************************************************************
long ignbin_mt(long n,double pp)
GENerate BINomial random deviate
Function
Generates a single random deviate from a binomial
distribution whose number of trials is N and whose
probability of an event in each trial is P.
Arguments
n  --> The number of trials in the binomial distribution
from which a random deviate is to be generated.
JJV (N >= 0)
pp --> The probability of an event in each trial of the
binomial distribution from which a random deviate
is to be generated.
JJV (0.0 <= PP <= 1.0)
ignbin <-- A random deviate yielding the number of events
from N independent trials, each of which has
a probability of event P.
Method
This is algorithm BTPE from:
Kachitvichyanukul, V. and Schmeiser, B. W.
Binomial Random Variate Generation.
Communications of the ACM, 31, 2
(February, 1988) 216.
**********************************************************************
SUBROUTINE BTPEC(N,PP,ISEED,JX)
BINOMIAL RANDOM VARIATE GENERATOR
MEAN .LT. 30 -- INVERSE CDF
MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
(SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
USABLE ALGORITHM.
REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
"BINOMIAL RANDOM VARIATE GENERATION,"
COMMUNICATIONS OF THE ACM, FORTHCOMING
WRITTEN:  SEPTEMBER 1980.
LAST REVISED:  MAY 1985, JULY 1987
REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
GENERATOR
ARGUMENTS
N : NUMBER OF BERNOULLI TRIALS            (INPUT)
PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
VARIABLES
PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
M:  INTEGER VALUE OF THE CURRENT MODE
FM:  FLOATING POINT VALUE OF THE CURRENT MODE
XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
P1:  AREA OF THE TRIANGLE
C:  HEIGHT OF THE PARALLELOGRAMS
XM:  CENTER OF THE TRIANGLE
XL:  LEFT END OF THE TRIANGLE
XR:  RIGHT END OF THE TRIANGLE
AL:  TEMPORARY VARIABLE
XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
P2:  AREA OF THE PARALLELOGRAMS
P3:  AREA OF THE LEFT EXPONENTIAL TAIL
P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
FROM THE REGION
V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
(REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
REJECT THE CANDIDATE VALUE
IX:  INTEGER CANDIDATE VALUE
X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
K:  ABSOLUTE VALUE OF (IX-M)
F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
ALSO USED IN THE INVERSE TRANSFORMATION
R: THE RATIO P/Q
G: CONSTANT USED IN CALCULATION OF PROBABILITY
MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
OF F WHEN IX IS GREATER THAN M
IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
CALCULATION OF F WHEN IX IS LESS THAN M
I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
YNORM: LOGARITHM OF NORMAL BOUND
ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
USED IN THE FINAL ACCEPT/REJECT TEST
QN: PROBABILITY OF NO SUCCESS IN N TRIALS
REMARK
IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
ARE NOT INVOLVED.
ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/
	{
	/* JJV changed initial values to ridiculous values */
	double psave = -1.0E37;
	long nsave = -214748365;
	long ignbin_mt,i,ix,ix1,k,m,mp,T1;
	double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
		x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

S10:
	/*
	*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
	JJV added checks to ensure 0.0 <= PP <= 1.0
	*/
	if(pp < 0.0) ERR_CRITICAL("PP < 0.0 in IGNBIN");
	if(pp > 1.0) ERR_CRITICAL("PP > 1.0 in IGNBIN");
	psave = pp;
	p = minF(psave,1.0-psave);
	q = 1.0-p;
S20:
	/*
	JJV added check to ensure N >= 0
	*/
	if(n < 0) ERR_CRITICAL("N < 0 in IGNBIN");
	xnp = n*p;
	nsave = n;
	if(xnp < 30.0) goto S140;
	ffm = xnp+p;
	m = ffm;
	fm = m;
	xnpq = xnp*q;
	p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
	xm = fm+0.5;
	xl = xm-p1;
	xr = xm+p1;
	c = 0.134+20.5/(15.3+fm);
	al = (ffm-xl)/(ffm-xl*p);
	xll = al*(1.0+0.5*al);
	al = (xr-ffm)/(xr*q);
	xlr = al*(1.0+0.5*al);
	p2 = p1*(1.0+c+c);
	p3 = p2+c/xll;
	p4 = p3+c/xlr;
S30:
	/*
	*****GENERATE VARIATE
	*/
	u = ranf_mt(state, tn)*p4;
	v = ranf_mt(state, tn);
	/*
	TRIANGULAR REGION
	*/
	if(u > p1) goto S40;
	ix = xm-p1*v+u;
	goto S170;
S40:
	/*
	PARALLELOGRAM REGION
	*/
	if(u > p2) goto S50;
	x = xl+(u-p1)/c;
	v = v*c+1.0 - ABS(xm-x)/p1;
	if(v > 1.0 || v <= 0.0) goto S30;
	ix = x;
	goto S70;
S50:
	/*
	LEFT TAIL
	*/
	if(u > p3) goto S60;
	ix = xl + log(v)/xll;
	if(ix < 0) goto S30;
	v *= ((u-p2)*xll);
	goto S70;
S60:
	/*
	RIGHT TAIL
	*/
	ix = xr - log(v)/xlr;
	if(ix > n) goto S30;
	v *= ((u-p3)*xlr);
S70:
	/*
	*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
	*/
	k = ABS(ix-m);
	if(k > 20 && k < xnpq/2-1) goto S130;
	/*
	EXPLICIT EVALUATION
	*/
	f = 1.0;
	r = p/q;
	g = (n+1)*r;
	T1 = m-ix;
	if(T1 < 0) goto S80;
	else if(T1 == 0) goto S120;
	else  goto S100;
S80:
	mp = m+1;
	for(i=mp; i<=ix; i++) f *= (g/i-r);
	goto S120;
S100:
	ix1 = ix+1;
	for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
	if(v <= f) goto S170;
	goto S30;
S130:
	/*
	SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
	*/
	amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
	ynorm = -(k*k/(2.0*xnpq));
	alv = log(v);
	if(alv < ynorm-amaxp) goto S170;
	if(alv > ynorm+amaxp) goto S30;
	/*
	STIRLING'S FORMULA TO MACHINE ACCURACY FOR
	THE FINAL ACCEPTANCE/REJECTION TEST
	*/
	x1 = ix+1.0;
	f1 = fm+1.0;
	z = n+1.0-fm;
	w = n-ix+1.0;
	z2 = z*z;
	x2 = x1*x1;
	f2 = f1*f1;
	w2 = w*w;
	if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
		(462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
		(132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
		(99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
		-140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
	goto S30;
S140:
	/*
	INVERSE CDF LOGIC FOR MEAN LESS THAN 30
	*/
	qn = pow(q,(double)n);
	r = p/q;
	g = r*(n+1);
S150:
	ix = 0;
	f = qn;
	u = ranf_mt(state, tn);
S160:
	if(u < f) goto S170;
	if(ix > 110) goto S150;
	u -= f;
	ix += 1;
	f *= (g/ix-r);
	goto S160;
S170:
	if(psave > 0.5) ix = n-ix;
	ignbin_mt = ix;
	return ignbin_mt;
}


// long ignpoi(double mu)
// /*
// **********************************************************************
//      long ignpoi(double mu)
//                     GENerate POIsson random deviate
//                               Function
//      Generates a single random deviate from a Poisson
//      distribution with mean MU.
//                               Arguments
//      mu --> The mean of the Poisson distribution from which
//             a random deviate is to be generated.
// 	    (mu >= 0.0)
//      ignpoi <-- The random deviate.
//                               Method
//      Renames KPOIS from TOMS as slightly modified by BWB to use RANF
//      instead of SUNIF.
//      For details see:
//                Ahrens, J.H. and Dieter, U.
//                Computer Generation of Poisson Deviates
//                From Modified Normal Distributions.
//                ACM Trans. Math. Software, 8, 2
//                (June 1982),163-179
// **********************************************************************
// **********************************************************************
                                                                      
                                                                      
//      P O I S S O N  DISTRIBUTION                                      
                                                                      
                                                                      
// **********************************************************************
// **********************************************************************
                                                                      
//      FOR DETAILS SEE:                                                 
                                                                      
//                AHRENS, J.H. AND DIETER, U.                            
//                COMPUTER GENERATION OF POISSON DEVIATES                
//                FROM MODIFIED NORMAL DISTRIBUTIONS.                    
//                ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. 
                                                                      
//      (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  
                                                                      
// **********************************************************************
//       INTEGER FUNCTION IGNPOI(IR,MU)
//      INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
//              MU=MEAN MU OF THE POISSON DISTRIBUTION
//      OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
//      MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
//      TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
//      COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
//      SEPARATION OF CASES A AND B
// */
// {
// extern double fsign( double num, double sign );
// static double a0 = -0.5;
// static double a1 = 0.3333333;
// static double a2 = -0.2500068;
// static double a3 = 0.2000118;
// static double a4 = -0.1661269;
// static double a5 = 0.1421878;
// static double a6 = -0.1384794;
// static double a7 = 0.125006;
// /* JJV changed the initial values of MUPREV and MUOLD */
// static double fact[10] = {
// 	1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
// 	};
// /* JJV added ll to the list, for Case A */
// long ignpoi,j,k,kflag,l,ll,m;
// double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,t,u,v,x,xx,pp[35];

//     if(mu < 10.0) goto S120;
// /*
//      C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
//      JJV changed l in Case A to ll
// */
//     s = sqrt(mu);
//     d = 6.0*mu*mu;
// /*
//              THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
//              PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
//              IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
// */
//     ll = (long) (mu-1.1484);
// S10:
// /*
//      STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
// */
//     g = mu+s*snorm();
//     if(g < 0.0) goto S20;
//     ignpoi = (long) (g);
// /*
//      STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
// */
//     if(ignpoi >= ll) return ignpoi;
// /*
//      STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
// */
//     fk = (double)ignpoi;
//     difmuk = mu-fk;
//     u = ranf();
//     if(d*u >= difmuk*difmuk*difmuk) return ignpoi;
// S20:
// /*
//      STEP P. PREPARATIONS FOR STEPS Q AND H.
//              (RECALCULATIONS OF PARAMETERS IF NECESSARY)
//              .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
//              THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
//              APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
//              C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
// */
//     omega = 0.3989423/s;
//     b1 = 4.166667E-2/mu;
//     b2 = 0.3*b1*b1;
//     c3 = 0.1428571*b1*b2;
//     c2 = b2-15.0*c3;
//     c1 = b1-6.0*b2+45.0*c3;
//     c0 = 1.0-b1+3.0*b2-15.0*c3;
//     c = 0.1069/mu;
// S30:
//     if(g < 0.0) goto S50;
// /*
//              'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
// */
//     kflag = 0;
//     goto S70;
// S40:
// /*
//      STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
// */
//     if(fy-u*fy <= py*exp(px-fx)) return ignpoi;
// S50:
// /*
//      STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
//              DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
//              (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
// */
//     e = sexpo();
//     u = ranf();
//     u += (u-1.0);
//     t = 1.8+fsign(e,u);
//     if(t <= -0.6744) goto S50;
//     ignpoi = (long) (mu+s*t);
//     fk = (double)ignpoi;
//     difmuk = mu-fk;
// /*
//              'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
// */
//     kflag = 1;
//     goto S70;
// S60:
// /*
//      STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
// */
//     if(c*fabs(u) > py*exp(px+e)-fy*exp(fx+e)) goto S50;
//     return ignpoi;
// S70:
// /*
//      STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
//              CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
// */
//     if(ignpoi >= 10) goto S80;
//     px = -mu;
//     py = pow(mu,(double)ignpoi)/ *(fact+ignpoi);
//     goto S110;
// S80:
// /*
//              CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
//              A0-A7 FOR ACCURACY WHEN ADVISABLE
//              .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
// */
//     del = 8.333333E-2/fk;
//     del -= (4.8*del*del*del);
//     v = difmuk/fk;
//     if(fabs(v) <= 0.25) goto S90;
//     px = fk*log(1.0+v)-difmuk-del;
//     goto S100;
// S90:
//     px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
// S100:
//     py = 0.3989423/sqrt(fk);
// S110:
//     x = (0.5-difmuk)/s;
//     xx = x*x;
//     fx = -0.5*xx;
//     fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
//     if(kflag <= 0) goto S40;
//     goto S60;
// S120:
// /*
//      C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
//      JJV changed MUPREV assignment to initial value
// */
//     m = maxF(1L,(long) (mu));
//     l = 0;
//     p = exp(-mu);
//     q = p0 = p;
// S130:
// /*
//      STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
// */
//     u = ranf();
//     ignpoi = 0;
//     if(u <= p0) return ignpoi;
// /*
//      STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
//              PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
//              (0.458=PP(9) FOR MU=10)
// */
//     if(l == 0) goto S150;
//     j = 1;
//     if(u > 0.458) j = minF(l,m);
//     for(k=j; k<=l; k++) {
//         if(u <= *(pp+k-1)) goto S180;
//     }
//     if(l == 35) goto S130;
// S150:
// /*
//      STEP C. CREATION OF NEW POISSON PROBABILITIES P
//              AND THEIR CUMULATIVES Q=PP(K)
// */
//     l += 1;
//     for(k=l; k<=35; k++) {
//         p = p*mu/(double)k;
//         q += p;
//         *(pp+k-1) = q;
//         if(u <= q) goto S170;
//     }
//     l = 35;
//     goto S130;
// S170:
//     l = k;
// S180:
//     ignpoi = k;
//     return ignpoi;
// }

// long ignpoi_mt(double mu, int tn)
// /*
// **********************************************************************
// long ignpoi_mt(double mu)
// GENerate POIsson random deviate
// Function
// Generates a single random deviate from a Poisson
// distribution with mean MU.
// Arguments
// mu --> The mean of the Poisson distribution from which
// a random deviate is to be generated.
// (mu >= 0.0)
// ignpoi_mt <-- The random deviate.
// Method
// Renames KPOIS from TOMS as slightly modified by BWB to use RANF
// instead of SUNIF.
// For details see:
// Ahrens, J.H. and Dieter, U.
// Computer Generation of Poisson Deviates
// From Modified Normal Distributions.
// ACM Trans. Math. Software, 8, 2
// (June 1982),163-179
// **********************************************************************
// **********************************************************************


// P O I S S O N  DISTRIBUTION                                      


// **********************************************************************
// **********************************************************************

// FOR DETAILS SEE:                                                 

// AHRENS, J.H. AND DIETER, U.                            
// COMPUTER GENERATION OF POISSON DEVIATES                
// FROM MODIFIED NORMAL DISTRIBUTIONS.                    
// ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. 

// (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  

// **********************************************************************
// INTEGER FUNCTION IGNPOI(IR,MU)
// INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
// MU=MEAN MU OF THE POISSON DISTRIBUTION
// OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
// MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
// TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
// COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
// SEPARATION OF CASES A AND B
// */
// 	{
// 	extern double fsign( double num, double sign );
// 	double a0 = -0.5;
// 	double a1 = 0.3333333;
// 	double a2 = -0.2500068;
// 	double a3 = 0.2000118;
// 	double a4 = -0.1661269;
// 	double a5 = 0.1421878;
// 	double a6 = -0.1384794;
// 	double a7 = 0.125006;
// 	/* JJV changed the initial values of MUPREV and MUOLD */
// 	double fact[10] = {
// 		1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
// 		};
// 	/* JJV added ll to the list, for Case A */
// 	long ignpoi_mt,j,k,kflag,l,ll,m;
// 	double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,t,u,v,x,xx,pp[35];

// 	if(mu < 10.0) goto S120;
// 	/*
// 	C A S E  A. (RECALCULATION OF S,D,LL IF MU HAS CHANGED)
// 	JJV changed l in Case A to ll
// 	*/
// 	s = sqrt(mu);
// 	d = 6.0*mu*mu;
// 	/*
// 	THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
// 	PROBABILITIES FK WHENEVER K >= M(MU). LL=IFIX(MU-1.1484)
// 	IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
// 	*/
// 	ll = (long) (mu-1.1484);
// S10:
// 	/*
// 	STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
// 	*/
// 	g = mu+s*snorm_mt(tn);
// 	if(g < 0.0) goto S20;
// 	ignpoi_mt = (long) (g);
// 	/*
// 	STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
// 	*/
// 	if(ignpoi_mt >= ll) return ignpoi_mt;
// 	/*
// 	STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
// 	*/
// 	fk = (double)ignpoi_mt;
// 	difmuk = mu-fk;
// 	u = ranf_mt(tn);
// 	if(d*u >= difmuk*difmuk*difmuk) return ignpoi_mt;
// S20:
// 	/*
// 	STEP P. PREPARATIONS FOR STEPS Q AND H.
// 	(RECALCULATIONS OF PARAMETERS IF NECESSARY)
// 	.3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
// 	THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
// 	APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
// 	C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
// 	*/
// 	omega = 0.3989423/s;
// 	b1 = 4.166667E-2/mu;
// 	b2 = 0.3*b1*b1;
// 	c3 = 0.1428571*b1*b2;
// 	c2 = b2-15.0*c3;
// 	c1 = b1-6.0*b2+45.0*c3;
// 	c0 = 1.0-b1+3.0*b2-15.0*c3;
// 	c = 0.1069/mu;
// S30:
// 	if(g < 0.0) goto S50;
// 	/*
// 	'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
// 	*/
// 	kflag = 0;
// 	goto S70;
// S40:
// 	/*
// 	STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
// 	*/
// 	if(fy-u*fy <= py*exp(px-fx)) return ignpoi_mt;
// S50:
// 	/*
// 	STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
// 	DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
// 	(IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
// 	*/
// 	e = sexpo_mt(tn);
// 	u = ranf_mt(tn);
// 	u += (u-1.0);
// 	t = 1.8+fsign(e,u);
// 	if(t <= -0.6744) goto S50;
// 	ignpoi_mt = (long) (mu+s*t);
// 	fk = (double)ignpoi_mt;
// 	difmuk = mu-fk;
// 	/*
// 	'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
// 	*/
// 	kflag = 1;
// 	goto S70;
// S60:
// 	/*
// 	STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
// 	*/
// 	if(c*fabs(u) > py*exp(px+e)-fy*exp(fx+e)) goto S50;
// 	return ignpoi_mt;
// S70:
// 	/*
// 	STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
// 	CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
// 	*/
// 	if(ignpoi_mt >= 10) goto S80;
// 	px = -mu;
// 	py = pow(mu,(double)ignpoi_mt)/ *(fact+ignpoi_mt);
// 	goto S110;
// S80:
// 	/*
// 	CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
// 	A0-A7 FOR ACCURACY WHEN ADVISABLE
// 	.8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
// 	*/
// 	del = 8.333333E-2/fk;
// 	del -= (4.8*del*del*del);
// 	v = difmuk/fk;
// 	if(fabs(v) <= 0.25) goto S90;
// 	px = fk*log(1.0+v)-difmuk-del;
// 	goto S100;
// S90:
// 	px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
// S100:
// 	py = 0.3989423/sqrt(fk);
// S110:
// 	x = (0.5-difmuk)/s;
// 	xx = x*x;
// 	fx = -0.5*xx;
// 	fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
// 	if(kflag <= 0) goto S40;
// 	goto S60;
// S120:
// 	/*
// 	C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
// 	JJV changed MUPREV assignment to initial value
// 	*/
// 	m = maxF(1L,(long) (mu));
// 	l = 0;
// 	p = exp(-mu);
// 	q = p0 = p;
// S130:
// 	/*
// 	STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
// 	*/
// 	u = ranf_mt(tn);
// 	ignpoi_mt = 0;
// 	if(u <= p0) return ignpoi_mt;
// 	/*
// 	STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
// 	PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
// 	(0.458=PP(9) FOR MU=10)
// 	*/
// 	if(l == 0) goto S150;
// 	j = 1;
// 	if(u > 0.458) j = minF(l,m);
// 	for(k=j; k<=l; k++) {
// 		if(u <= *(pp+k-1)) goto S180;
// 		}
// 	if(l == 35) goto S130;
// S150:
// /*
// STEP C. CREATION OF NEW POISSON PROBABILITIES P
// AND THEIR CUMULATIVES Q=PP(K)
// */
// 	l += 1;
// 	for(k=l; k<=35; k++) {
// 		p = p*mu/(double)k;
// 		q += p;
// 		*(pp+k-1) = q;
// 		if(u <= q) goto S170;
// 		}
// 	l = 35;
// 	goto S130;
// 	S170:
// 	l = k;
// 	S180:
// 	ignpoi_mt = k;
// 	return ignpoi_mt;
// }


// double sexpo(void)
// /*
// **********************************************************************
                                                                      
                                                                      
//      (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                
                                                                      
                                                                      
// **********************************************************************
// **********************************************************************
                                                                      
//      FOR DETAILS SEE:                                                 
                                                                      
//                AHRENS, J.H. AND DIETER, U.                            
//                COMPUTER METHODS FOR SAMPLING FROM THE                 
//                EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
//                COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               
                                                                      
//      ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
//      'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
                                                                      
//      Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
//      SUNIF.  The argument IR thus goes away.                          
                                                                      
// **********************************************************************
//      Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
//      (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
// */
// {
// static double q[8] = {
//     0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
// 	0.99999999999999989
// };
// long i;
// double sexpo,a,u,ustar,umin;

// 	a = 0.0;
//     u = ranf();
//     goto S30;
// S20:
//     a += q[0];
// S30:
//     u += u;
// /*
//  * JJV changed the following to reflect the true algorithm and prevent
//  * JJV unpredictable behavior if U is initially 0.5.
//  *  if(u <= 1.0) goto S20;
//  */
//     if(u < 1.0) goto S20;
//     u -= 1.0;
//     if(u > q[0]) goto S60;
//     sexpo = a+u;
//     return sexpo;
// S60:
//     i = 1;
//     ustar = ranf();
//     umin = ustar;
// S70:
//     ustar = ranf();
//     if(ustar < umin) umin = ustar;
//     i += 1;
//     if(u > q[i-1]) goto S70;
//     return  a+umin*q[0];
// }

// double sexpo_mt(int tn)
// /*
// **********************************************************************


// (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                


// **********************************************************************
// **********************************************************************

// FOR DETAILS SEE:                                                 

// AHRENS, J.H. AND DIETER, U.                            
// COMPUTER METHODS FOR SAMPLING FROM THE                 
// EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
// COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               

// ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
// 'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       

// Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
// SUNIF.  The argument IR thus goes away.                          

// **********************************************************************
// Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
// (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
// */
// 	{
// 	double q[8] = {
// 		0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,
// 			0.99999999999999989
// 		};
// 	long i;
// 	double sexpo_mt,a,u,ustar,umin;

// 	a = 0.0;
// 	u = ranf_mt(tn);
// 	goto S30;
// S20:
// 	a += q[0];
// S30:
// 	u += u;
// 	/*
// 	* JJV changed the following to reflect the true algorithm and prevent
// 	* JJV unpredictable behavior if U is initially 0.5.
// 	*  if(u <= 1.0) goto S20;
// 	*/
// 	if(u < 1.0) goto S20;
// 	u -= 1.0;
// 	if(u > q[0]) goto S60;
// 	sexpo_mt = a+u;
// 	return sexpo_mt;
// S60:
// 	i = 1;
// 	ustar = ranf_mt(tn);
// 	umin = ustar;
// S70:
// 	ustar = ranf_mt(tn);
// 	if(ustar < umin) umin = ustar;
// 	i += 1;
// 	if(u > q[i-1]) goto S70;
// 	return  a+umin*q[0];
// }

// double snorm(void)
// /*
// **********************************************************************
                                                                      
                                                                      
//      (STANDARD-)  N O R M A L  DISTRIBUTION                           
                                                                      
                                                                      
// **********************************************************************
// **********************************************************************
                                                                      
//      FOR DETAILS SEE:                                                 
                                                                      
//                AHRENS, J.H. AND DIETER, U.                            
//                EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
//                SAMPLING FROM THE NORMAL DISTRIBUTION.                 
//                MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          
                                                                      
//      ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
//      (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  
                                                                      
//      Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
//      SUNIF.  The argument IR thus goes away.                          
                                                                      
// **********************************************************************
//      THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
//      H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
// */
// {
// static double a[32] = {
//     0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
//     0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
//     0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
//     1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
//     1.862732,2.153875
// };
// static double d[31] = {
//     0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
//     0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
//     0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
//     0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
// };
// static double t[31] = {
//     7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
//     1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
//     2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
//     4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
//     9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
// };
// static double h[31] = {
//     3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
//     4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
//     4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
//     5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
//     8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
// };
// 	long i; //made this non-static: ggilani 27/11/14
// 	double snorm,u,s,ustar,aa,w,y,tt; //made this non-static: ggilani 27/11/14
//     u = ranf();
//     s = 0.0;
//     if(u > 0.5) s = 1.0;
//     u += (u-s);
//     u = 32.0*u;
//     i = (long) (u);
//     if(i == 32) i = 31;
//     if(i == 0) goto S100;
// /*
//                                 START CENTER
// */
//     ustar = u-(double)i;
//     aa = *(a+i-1);
// S40:
//     if(ustar <= *(t+i-1)) goto S60;
//     w = (ustar-*(t+i-1))**(h+i-1);
// S50:
// /*
//                                 EXIT   (BOTH CASES)
// */
//     y = aa+w;
//     snorm = y;
//     if(s == 1.0) snorm = -y;
//     return snorm;
// S60:
// /*
//                                 CENTER CONTINUED
// */
//     u = ranf();
//     w = u*(*(a+i)-aa);
//     tt = (0.5*w+aa)*w;
//     goto S80;
// S70:
//     tt = u;
//     ustar = ranf();
// S80:
//     if(ustar > tt) goto S50;
//     u = ranf();
//     if(ustar >= u) goto S70;
//     ustar = ranf();
//     goto S40;
// S100:
// /*
//                                 START TAIL
// */
//     i = 6;
//     aa = *(a+31);
//     goto S120;
// S110:
//     aa += *(d+i-1);
//     i += 1;
// S120:
//     u += u;
//     if(u < 1.0) goto S110;
//     u -= 1.0;
// S140:
//     w = u**(d+i-1);
//     tt = (0.5*w+aa)*w;
//     goto S160;
// S150:
//     tt = u;
// S160:
//     ustar = ranf();
//     if(ustar > tt) goto S50;
//     u = ranf();
//     if(ustar >= u) goto S150;
//     u = ranf();
//     goto S140;
// }

// double snorm_mt(int tn)
// /*
// **********************************************************************


// (STANDARD-)  N O R M A L  DISTRIBUTION                           


// **********************************************************************
// **********************************************************************

// FOR DETAILS SEE:                                                 

// AHRENS, J.H. AND DIETER, U.                            
// EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
// SAMPLING FROM THE NORMAL DISTRIBUTION.                 
// MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          

// ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
// (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  

// Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
// SUNIF.  The argument IR thus goes away.                          

// **********************************************************************
// THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
// H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
// */
// 	{
// 	static double a[32] = {
// 		0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
// 			0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
// 			0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
// 			1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
// 			1.862732,2.153875
// 		};
// 	static double d[31] = {
// 		0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
// 			0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
// 			0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
// 			0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
// 		};
// 	static double t[31] = {
// 		7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
// 			1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
// 			2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
// 			4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
// 			9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
// 		};
// 	static double h[31] = {
// 		3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
// 			4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
// 			4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
// 			5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
// 			8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
// 		};
// 	long i;
// 	double snorm_mt,u,s,ustar,aa,w,y,tt;
// 	u = ranf_mt(tn);
// 	s = 0.0;
// 	if(u > 0.5) s = 1.0;
// 	u += (u-s);
// 	u = 32.0*u;
// 	i = (long) (u);
// 	if(i == 32) i = 31;
// 	if(i == 0) goto S100;
// 	/*
// 	START CENTER
// 	*/
// 	ustar = u-(double)i;
// 	aa = *(a+i-1);
// S40:
// 	if(ustar <= *(t+i-1)) goto S60;
// 	w = (ustar-*(t+i-1))**(h+i-1);
// S50:
// 	/*
// 	EXIT   (BOTH CASES)
// 	*/
// 	y = aa+w;
// 	snorm_mt = y;
// 	if(s == 1.0) snorm_mt = -y;
// 	return snorm_mt;
// S60:
// 	/*
// 	CENTER CONTINUED
// 	*/
// 	u = ranf_mt(tn);
// 	w = u*(*(a+i)-aa);
// 	tt = (0.5*w+aa)*w;
// 	goto S80;
// S70:
// 	tt = u;
// 	ustar = ranf_mt(tn);
// S80:
// 	if(ustar > tt) goto S50;
// 	u = ranf_mt(tn);
// 	if(ustar >= u) goto S70;
// 	ustar = ranf_mt(tn);
// 	goto S40;
// S100:
// 	/*
// 	START TAIL
// 	*/
// 	i = 6;
// 	aa = *(a+31);
// 	goto S120;
// S110:
// 	aa += *(d+i-1);
// 	i += 1;
// S120:
// 	u += u;
// 	if(u < 1.0) goto S110;
// 	u -= 1.0;
// S140:
// 	w = u**(d+i-1);
// 	tt = (0.5*w+aa)*w;
// 	goto S160;
// S150:
// 	tt = u;
// S160:
// 	ustar = ranf_mt(tn);
// 	if(ustar > tt) goto S50;
// 	u = ranf_mt(tn);
// 	if(ustar >= u) goto S150;
// 	u = ranf_mt(tn);
// 	goto S140;
// }
// double fsign( double num, double sign )
// /* Transfers sign of argument sign to argument num */
// {
// if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
//     return -num;
// else return num;
// }

// double gengam(double a,double r)
// /*
// **********************************************************************
// double gengam(double a,double r)
// GENerates random deviates from GAMma distribution
// Function
// Generates random deviates from the gamma distribution whose
// density is
// (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
// Arguments
// a --> Location parameter of Gamma distribution
// JJV   (a > 0)
// r --> Shape parameter of Gamma distribution
// JJV   (r > 0)
// Method
// Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
// instead of SUNIF.
// For details see:
// (Case R >= 1.0)
// Ahrens, J.H. and Dieter, U.
// Generating Gamma Variates by a
// Modified Rejection Technique.
// Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
// Algorithm GD
// JJV altered following to reflect argument ranges
// (Case 0.0 < R < 1.0)
// Ahrens, J.H. and Dieter, U.
// Computer Methods for Sampling from Gamma,
// Beta, Poisson and Binomial Distributions.
// Computing, 12 (1974), 223-246/
// Adapted algorithm GS.
// **********************************************************************
// */
// {
// 	static double gengam;
// 	/* JJV added argument checker */
// 	if(a > 0.0 && r > 0.0) goto S10;
// 	fputs(" A or R nonpositive in GENGAM - abort!\n",stderr);
// 	fprintf(stderr," A value: %16.6E R value: %16.6E\n",a,r);
// 	exit(1);
// S10:
// 	gengam = sgamma(r);
// 	gengam /= a;
// 	return gengam;
// }

// double sgamma(double a)
// /*
// **********************************************************************


// (STANDARD-)  G A M M A  DISTRIBUTION                             


// **********************************************************************
// **********************************************************************

// PARAMETER  A >= 1.0  !                                 

// **********************************************************************

// FOR DETAILS SEE:                                                 

// AHRENS, J.H. AND DIETER, U.                            
// GENERATING GAMMA VARIATES BY A                         
// MODIFIED REJECTION TECHNIQUE.                          
// COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  

// STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
// (STRAIGHTFORWARD IMPLEMENTATION)     

// Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
// SUNIF.  The argument IR thus goes away.                          

// **********************************************************************

// PARAMETER  0.0 < A < 1.0  !                            

// **********************************************************************

// FOR DETAILS SEE:                                                 

// AHRENS, J.H. AND DIETER, U.                            
// COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
// BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
// COMPUTING, 12 (1974), 223 - 246.                       

// (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    

// **********************************************************************
// INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
// OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
// COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
// COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
// COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
// PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
// SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
// */
// {
// 	extern double fsign( double num, double sign );
// 	static double q1 = 4.166669E-2;
// 	static double q2 = 2.083148E-2;
// 	static double q3 = 8.01191E-3;
// 	static double q4 = 1.44121E-3;
// 	static double q5 = -7.388E-5;
// 	static double q6 = 2.4511E-4;
// 	static double q7 = 2.424E-4;
// 	static double a1 = 0.3333333;
// 	static double a2 = -0.250003;
// 	static double a3 = 0.2000062;
// 	static double a4 = -0.1662921;
// 	static double a5 = 0.1423657;
// 	static double a6 = -0.1367177;
// 	static double a7 = 0.1233795;
// 	static double e1 = 1.0;
// 	static double e2 = 0.4999897;
// 	static double e3 = 0.166829;
// 	static double e4 = 4.07753E-2;
// 	static double e5 = 1.0293E-2;
// 	static double aa = 0.0;
// 	static double aaa = 0.0;
// 	static double sqrt32 = 5.656854;
// 	/* JJV added b0 to fix rare and subtle bug */
// 	static double sgamma,s2,s,d,t,x,u,r,q0,b,b0,si,c,v,q,e,w,p; //remove static from these variables?
// 	if(a == aa) goto S10;
// 	if(a < 1.0) goto S120;
// 	/*
// 	STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
// 	*/
// 	aa = a;
// 	s2 = a-0.5;
// 	s = sqrt(s2);
// 	d = sqrt32-12.0*s;
// S10:
// 	/*
// 	STEP  2:  T=STANDARD NORMAL DEVIATE,
// 	X=(S,1/2)-NORMAL DEVIATE.
// 	IMMEDIATE ACCEPTANCE (I)
// 	*/
// 	t = snorm();
// 	x = s+0.5*t;
// 	sgamma = x*x;
// 	if(t >= 0.0) return sgamma;
// 	/*
// 	STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
// 	*/
// 	u = ranf();
// 	if(d*u <= t*t*t) return sgamma;
// 	/*
// 	STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
// 	*/
// 	if(a == aaa) goto S40;
// 	aaa = a;
// 	r = 1.0/ a;
// 	q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
// 	/*
// 	APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
// 	THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
// 	C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
// 	*/
// 	if(a <= 3.686) goto S30;
// 	if(a <= 13.022) goto S20;
// 	/*
// 	CASE 3:  A .GT. 13.022
// 	*/
// 	b = 1.77;
// 	si = 0.75;
// 	c = 0.1515/s;
// 	goto S40;
// S20:
// 	/*
// 	CASE 2:  3.686 .LT. A .LE. 13.022
// 	*/
// 	b = 1.654+7.6E-3*s2;
// 	si = 1.68/s+0.275;
// 	c = 6.2E-2/s+2.4E-2;
// 	goto S40;
// S30:
// 	/*
// 	CASE 1:  A .LE. 3.686
// 	*/
// 	b = 0.463+s+0.178*s2;
// 	si = 1.235;
// 	c = 0.195/s-7.9E-2+1.6E-1*s;
// S40:
// 	/*
// 	STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
// 	*/
// 	if(x <= 0.0) goto S70;
// 	/*
// 	STEP  6:  CALCULATION OF V AND QUOTIENT Q
// 	*/
// 	v = t/(s+s);
// 	if(fabs(v) <= 0.25) goto S50;
// 	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
// 	goto S60;
// S50:
// 	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
// S60:
// 	/*
// 	STEP  7:  QUOTIENT ACCEPTANCE (Q)
// 	*/
// 	if(log(1.0-u) <= q) return sgamma;
// S70:
// 	/*
// 	STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
// 	U= 0,1 -UNIFORM DEVIATE
// 	T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
// 	*/
// 	e = sexpo();
// 	u = ranf();
// 	u += (u-1.0);
// 	t = b+fsign(si*e,u);
// 	/*
// 	STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
// 	*/
// 	if(t < -0.7187449) goto S70;
// 	/*
// 	STEP 10:  CALCULATION OF V AND QUOTIENT Q
// 	*/
// 	v = t/(s+s);
// 	if(fabs(v) <= 0.25) goto S80;
// 	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
// 	goto S90;
// S80:
// 	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
// S90:
// 	/*
// 	STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
// 	*/
// 	if(q <= 0.0) goto S70;
// 	if(q <= 0.5) goto S100;
// 	/*
// 	* JJV modified the code through line 115 to handle large Q case
// 	*/
// 	if(q < 15.0) goto S95;
// 	/*
// 	* JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
// 	* JJV so reformulate test at 110 in terms of one EXP, if not too big
// 	* JJV 87.49823 is close to the largest real which can be
// 	* JJV exponentiated (87.49823 = log(1.0E38))
// 	*/
// 	if((q+e-0.5*t*t) > 87.49823) goto S115;
// 	if(c*fabs(u) > exp(q+e-0.5*t*t)) goto S70;
// 	goto S115;
// S95:
// 	w = exp(q)-1.0;
// 	goto S110;
// S100:
// 	w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
// S110:
// 	/*
// 	IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
// 	*/
// 	if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
// S115:
// 	x = s+0.5*t;
// 	sgamma = x*x;
// 	return sgamma;
// S120:
// 	/*
// 	ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))

// 	JJV changed B to B0 (which was added to declarations for this)
// 	JJV in 120 to END to fix rare and subtle bug.
// 	JJV Line: 'aa = 0.0' was removed (unnecessary, wasteful).
// 	JJV Reasons: the state of AA only serves to tell the A >= 1.0
// 	JJV case if certain A-dependent constants need to be recalculated.
// 	JJV The A < 1.0 case (here) no longer changes any of these, and
// 	JJV the recalculation of B (which used to change with an
// 	JJV A < 1.0 call) is governed by the state of AAA anyway.
// 	aa = 0.0;
// 	*/
// 	b0 = 1.0+0.3678794*a;
// S130:
// 	p = b0*ranf();
// 	if(p >= 1.0) goto S140;
// 	sgamma = exp(log(p)/ a);
// 	if(sexpo() < sgamma) goto S130;
// 	return sgamma;
// S140:
// 	sgamma = -log((b0-p)/ a);
// 	if(sexpo() < (1.0-a)*log(sgamma)) goto S130;
// 	return sgamma;
// }

// double gengam_mt(double a,double r,int tn)
// /*
// **********************************************************************
// double gengam(double a,double r)
// GENerates random deviates from GAMma distribution
// Function
// Generates random deviates from the gamma distribution whose
// density is
// (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
// Arguments
// a --> Location parameter of Gamma distribution
// JJV   (a > 0)
// r --> Shape parameter of Gamma distribution
// JJV   (r > 0)
// Method
// Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
// instead of SUNIF.
// For details see:
// (Case R >= 1.0)
// Ahrens, J.H. and Dieter, U.
// Generating Gamma Variates by a
// Modified Rejection Technique.
// Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
// Algorithm GD
// JJV altered following to reflect argument ranges
// (Case 0.0 < R < 1.0)
// Ahrens, J.H. and Dieter, U.
// Computer Methods for Sampling from Gamma,
// Beta, Poisson and Binomial Distributions.
// Computing, 12 (1974), 223-246/
// Adapted algorithm GS.
// **********************************************************************
// */
// 	{
// 	static double gengam_mt;
// 	/* JJV added argument checker */
// 	if(a > 0.0 && r > 0.0) goto S10;
// 	fputs(" A or R nonpositive in GENGAM - abort!\n",stderr);
// 	fprintf(stderr," A value: %16.6E R value: %16.6E\n",a,r);
// 	exit(1);
// S10:
// 	gengam_mt = sgamma_mt(r,tn);
// 	gengam_mt /= a;
// 	return gengam_mt;
// }

// double sgamma_mt(double a,int tn)
// /*
// **********************************************************************


// (STANDARD-)  G A M M A  DISTRIBUTION                             


// **********************************************************************
// **********************************************************************

// PARAMETER  A >= 1.0  !                                 

// **********************************************************************

// FOR DETAILS SEE:                                                 

// AHRENS, J.H. AND DIETER, U.                            
// GENERATING GAMMA VARIATES BY A                         
// MODIFIED REJECTION TECHNIQUE.                          
// COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  

// STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
// (STRAIGHTFORWARD IMPLEMENTATION)     

// Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
// SUNIF.  The argument IR thus goes away.                          

// **********************************************************************

// PARAMETER  0.0 < A < 1.0  !                            

// **********************************************************************

// FOR DETAILS SEE:                                                 

// AHRENS, J.H. AND DIETER, U.                            
// COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
// BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
// COMPUTING, 12 (1974), 223 - 246.                       

// (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    

// **********************************************************************
// INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
// OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
// COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
// COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
// COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
// PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
// SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
// */
// 	{
// 	extern double fsign( double num, double sign );
// 	static double q1 = 4.166669E-2;
// 	static double q2 = 2.083148E-2;
// 	static double q3 = 8.01191E-3;
// 	static double q4 = 1.44121E-3;
// 	static double q5 = -7.388E-5;
// 	static double q6 = 2.4511E-4;
// 	static double q7 = 2.424E-4;
// 	static double a1 = 0.3333333;
// 	static double a2 = -0.250003;
// 	static double a3 = 0.2000062;
// 	static double a4 = -0.1662921;
// 	static double a5 = 0.1423657;
// 	static double a6 = -0.1367177;
// 	static double a7 = 0.1233795;
// 	static double e1 = 1.0;
// 	static double e2 = 0.4999897;
// 	static double e3 = 0.166829;
// 	static double e4 = 4.07753E-2;
// 	static double e5 = 1.0293E-2;
// 	static double aa = 0.0;
// 	static double aaa = 0.0;
// 	static double sqrt32 = 5.656854;
// 	/* JJV added b0 to fix rare and subtle bug */
// 	static double sgamma_mt,s2,s,d,t,x,u,r,q0,b,b0,si,c,v,q,e,w,p; 
// 	if(a == aa) goto S10;
// 	if(a < 1.0) goto S120;
// 	/*
// 	STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
// 	*/
// 	aa = a;
// 	s2 = a-0.5;
// 	s = sqrt(s2);
// 	d = sqrt32-12.0*s;
// S10:
// 	/*
// 	STEP  2:  T=STANDARD NORMAL DEVIATE,
// 	X=(S,1/2)-NORMAL DEVIATE.
// 	IMMEDIATE ACCEPTANCE (I)
// 	*/
// 	t = snorm_mt(tn); //make this multi-threaded? 28/11/14
// 	x = s+0.5*t;
// 	sgamma_mt = x*x;
// 	if(t >= 0.0) return sgamma_mt;
// 	/*
// 	STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
// 	*/
// 	u = ranf_mt(tn);
// 	if(d*u <= t*t*t) return sgamma_mt;
// 	/*
// 	STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
// 	*/
// 	if(a == aaa) goto S40;
// 	aaa = a;
// 	r = 1.0/ a;
// 	q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
// 	/*
// 	APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
// 	THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
// 	C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
// 	*/
// 	if(a <= 3.686) goto S30;
// 	if(a <= 13.022) goto S20;
// 	/*
// 	CASE 3:  A .GT. 13.022
// 	*/
// 	b = 1.77;
// 	si = 0.75;
// 	c = 0.1515/s;
// 	goto S40;
// S20:
// 	/*
// 	CASE 2:  3.686 .LT. A .LE. 13.022
// 	*/
// 	b = 1.654+7.6E-3*s2;
// 	si = 1.68/s+0.275;
// 	c = 6.2E-2/s+2.4E-2;
// 	goto S40;
// S30:
// 	/*
// 	CASE 1:  A .LE. 3.686
// 	*/
// 	b = 0.463+s+0.178*s2;
// 	si = 1.235;
// 	c = 0.195/s-7.9E-2+1.6E-1*s;
// S40:
// 	/*
// 	STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
// 	*/
// 	if(x <= 0.0) goto S70;
// 	/*
// 	STEP  6:  CALCULATION OF V AND QUOTIENT Q
// 	*/
// 	v = t/(s+s);
// 	if(fabs(v) <= 0.25) goto S50;
// 	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
// 	goto S60;
// S50:
// 	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
// S60:
// 	/*
// 	STEP  7:  QUOTIENT ACCEPTANCE (Q)
// 	*/
// 	if(log(1.0-u) <= q) return sgamma_mt;
// S70:
// 	/*
// 	STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
// 	U= 0,1 -UNIFORM DEVIATE
// 	T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
// 	*/
// 	e = sexpo_mt(tn);
// 	u = ranf_mt(tn);
// 	u += (u-1.0);
// 	t = b+fsign(si*e,u);
// 	/*
// 	STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
// 	*/
// 	if(t < -0.7187449) goto S70;
// 	/*
// 	STEP 10:  CALCULATION OF V AND QUOTIENT Q
// 	*/
// 	v = t/(s+s);
// 	if(fabs(v) <= 0.25) goto S80;
// 	q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
// 	goto S90;
// S80:
// 	q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
// S90:
// 	/*
// 	STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
// 	*/
// 	if(q <= 0.0) goto S70;
// 	if(q <= 0.5) goto S100;
// 	/*
// 	* JJV modified the code through line 115 to handle large Q case
// 	*/
// 	if(q < 15.0) goto S95;
// 	/*
// 	* JJV Here Q is large enough that Q = log(exp(Q) - 1.0) (for real Q)
// 	* JJV so reformulate test at 110 in terms of one EXP, if not too big
// 	* JJV 87.49823 is close to the largest real which can be
// 	* JJV exponentiated (87.49823 = log(1.0E38))
// 	*/
// 	if((q+e-0.5*t*t) > 87.49823) goto S115;
// 	if(c*fabs(u) > exp(q+e-0.5*t*t)) goto S70;
// 	goto S115;
// S95:
// 	w = exp(q)-1.0;
// 	goto S110;
// S100:
// 	w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
// S110:
// 	/*
// 	IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
// 	*/
// 	if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
// S115:
// 	x = s+0.5*t;
// 	sgamma_mt = x*x;
// 	return sgamma_mt;
// S120:
// 	/*
// 	ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))

// 	JJV changed B to B0 (which was added to declarations for this)
// 	JJV in 120 to END to fix rare and subtle bug.
// 	JJV Line: 'aa = 0.0' was removed (unnecessary, wasteful).
// 	JJV Reasons: the state of AA only serves to tell the A >= 1.0
// 	JJV case if certain A-dependent constants need to be recalculated.
// 	JJV The A < 1.0 case (here) no longer changes any of these, and
// 	JJV the recalculation of B (which used to change with an
// 	JJV A < 1.0 call) is governed by the state of AAA anyway.
// 	aa = 0.0;
// 	*/
// 	b0 = 1.0+0.3678794*a;
// S130:
// 	p = b0*ranf_mt(tn);
// 	if(p >= 1.0) goto S140;
// 	sgamma_mt = exp(log(p)/ a);
// 	if(sexpo_mt(tn) < sgamma_mt) goto S130;
// 	return sgamma_mt;
// S140:
// 	sgamma_mt = -log((b0-p)/ a);
// 	if(sexpo_mt(tn) < (1.0-a)*log(sgamma_mt)) goto S130;
// 	return sgamma_mt;
// 	}


// /* function gen_snorm
//  * purpose: my own implementation of sampling from a uniform distribution, using Marsaglia polar method
//  *
//  * author: ggilani, date: 28/11/14
//  */
// double gen_norm(double mu, double sd)
// {
// 	double u,v,x,y,S;

// 	do
// 	{
// 		u=2*ranf()-1; //u and v are uniform random numbers on the interval [-1,1]
// 		v=2*ranf()-1;

// 		//calculate S=U^2+V^2
// 		S=u*u+v*v;
// 	}
// 	while(S>=1||S==0);

// 	//calculate x,y - both of which are normally distributed
// 	x=u*sqrt((-2*log(S))/S);
// 	y=v*sqrt((-2*log(S))/S);

// 	//return x
// 	return x*sd+mu;
// }

// /*function gen_snorm
//  * purpose: my own implementation of sampling from a uniform distribution, using Marsaglia polar method, but for multi-threading
//  *
//  * author: ggilani, date: 28/11/14
//  */
// double gen_norm_mt(double mu, double sd, int tn)
// {
// 	double u,v,x,y,S;

// 	do
// 	{
// 		u=2*ranf_mt(tn)-1; //u and v are uniform random numbers on the interval [-1,1]
// 		v=2*ranf_mt(tn)-1;

// 		//calculate S=U^2+V^2
// 		S=u*u+v*v;
// 	}
// 	while(S>=1||S==0);

// 	//calculate x,y - both of which are normally distributed
// 	x=u*sqrt((-2*log(S))/S);
// 	y=v*sqrt((-2*log(S))/S);

// 	//return x
// 	return x*sd+mu;
// }

// /*function gen_gamma
//  * purpose: my own implementation of sampling from a gamma distribution, using Marsaglia-Tsang method
//  *
//  * author: ggilani, date: 01/12/14
//  */
// double gen_gamma(double beta,double alpha)
// {
// 	double d,c,u,v,z,f,alpha2,gamma;

// 	//error statment if either beta or alpha are <=0, as gamma distribution is undefined in this case
// 	if((beta<=0)||(alpha<=0))
// 	{
// 		ERR_CRITICAL("Gamma distribution parameters in undefined range!\n");
// 	}

// 	//method is slightly different depending on whether alpha is greater than or equal to 1, or less than 1
// 	if(alpha>=1)
// 	{
// 		d=alpha-(1.0/3.0);
// 		c=1.0/(3.0*sqrt(d));
// 		do
// 		{
// 			//sample one random number from uniform distribution and one from standard normal distribution
// 			u=ranf();
// 			z=gen_norm(0,1);
// 			v=1+z*c;
// 			v=v*v*v;
// 		}
// 		while((z<=(-1.0/c))||(log(u)>=(0.5*z*z+d-d*v+d*log(v))));
// 		//if beta is equal to 1, there is no scale. If beta is not equal to one, return scaled value
// 		if(beta!=1)
// 		{
// 			return (d*v)/beta;
// 		}
// 		else
// 		{
// 			return d*v;
// 		}
// 	}
// 	else
// 	{
// 		//if alpha is less than 1, initially sample from gamma(beta,alpha+1)
// 		alpha2=alpha+1;
// 		d=alpha2-(1.0/3.0);
// 		c=1.0/(3.0*sqrt(d));
// 		do
// 		{
// 			//sample one random number from uniform distribution and one from standard normal distribution
// 			u=ranf();
// 			z=gen_norm(0,1);
// 			v=1+z*c;
// 			v=v*v*v;
// 		}
// 		while((z<=(-1.0/c))||(log(u)>=(0.5*z*z+d-d*v+d*log(v))));
// 		//if beta is equal to 1, there is no scale. If beta is not equal to one, return scaled value
// 		if(beta!=1)
// 		{
// 			gamma=(d*v)/beta;
// 		}
// 		else
// 		{
// 			gamma=d*v;
// 		}
// 		//now rescale again to take into account that alpha is less than 1
// 		f=pow(ranf(),(1.0/alpha));
// 		//return gamma scaled by f
// 		return gamma*f;
// 	}

// }

// /*function gen_gamma_mt
//  * purpose: my own implementation of sampling from a gamma distribution, using Marsaglia-Tsang method, but for multi-threading
//  *
//  * author: ggilani, date: 01/12/14
//  */
// double gen_gamma_mt(double beta,double alpha,int tn)
// {
// 	double d,c,u,v,z,f,alpha2,gamma;

// 	//error statment if either beta or alpha are <=0, as gamma distribution is undefined in this case
// 	if((beta<=0)||(alpha<=0))
// 	{
// 		ERR_CRITICAL("Gamma distribution parameters in undefined range!\n");
// 	}

// 	//method is slightly different depending on whether alpha is greater than or equal to 1, or less than 1
// 	if(alpha>=1)
// 	{
// 		d=alpha-(1.0/3.0);
// 		c=1.0/(3.0*sqrt(d));
// 		do
// 		{
// 			//sample one random number from uniform distribution and one from standard normal distribution
// 			u=ranf_mt(tn);
// 			z=gen_norm_mt(0,1,tn);
// 			v=1+z*c;
// 			v=v*v*v;
// 		}
// 		while((z<=(-1.0/c))||(log(u)>=(0.5*z*z+d-d*v+d*log(v))));
// 		//if beta is equal to 1, there is no scale. If beta is not equal to one, return scaled value
// 		if(beta!=1)
// 		{
// 			return (d*v)/beta;
// 		}
// 		else
// 		{
// 			return d*v;
// 		}
// 	}
// 	else
// 	{
// 		//if alpha is less than 1, initially sample from gamma(beta,alpha+1)
// 		alpha2=alpha+1;
// 		d=alpha2-(1.0/3.0);
// 		c=1.0/(3.0*sqrt(d));
// 		do
// 		{
// 			//sample one random number from uniform distribution and one from standard normal distribution
// 			u=ranf_mt(tn);
// 			z=gen_norm_mt(0,1,tn);
// 			v=1+z*c;
// 			v=v*v*v;
// 		}
// 		while((z<=(-1.0/c))||(log(u)>=(0.5*z*z+d-d*v+d*log(v))));
// 		//if beta is equal to 1, there is no scale. If beta is not equal to one, return scaled value
// 		if(beta!=1)
// 		{
// 			gamma=(d*v)/beta;
// 		}
// 		else
// 		{
// 			gamma=d*v;
// 		}
// 		//now rescale again to take into account that alpha is less than 1
// 		f=pow(ranf_mt(tn),(1.0/alpha));
// 		//return gamma scaled by f
// 		return gamma*f;
// 	}
// }


// /* function genbeta(double alpha,double beta)
//  * purpose: to generate samples from a beta distribution with parameters alpha and beta
//  *
//  * parameters:
//  *	shape parameter alpha
//  *	shape parameter beta
//  *
//  * returns: double from the beta distribution
//  *
//  * author: ggilani, date: 26/11/14
//  */
// double gen_beta(double alpha, double beta)
// {
// 	double randgam_a,randgam_b;

// 	//error statment if either beta or alpha are <=0, as gamma distribution is undefined in this case
// 	if((beta<=0)||(alpha<=0))
// 	{
// 		ERR_CRITICAL("Beta distribution parameters in undefined range!\n");
// 	}

// 	randgam_a=gen_gamma(1,alpha);
// 	randgam_b=gen_gamma(1,beta);

// 	return randgam_a/(randgam_a+randgam_b);
// }

// /* function genbeta_mt(double alpha,double beta, int tn)
//  * purpose: to generate multi thread samples from a beta distribution with parameters alpha and beta
//  *
//  * parameters:
//  *	shape parameter alpha
//  *	shape parameter beta
//  *	thread number tn
//  *
//  * returns: double from the beta distribution
//  *
//  * author: ggilani, date: 26/11/14
//  */
// double gen_beta_mt(double alpha, double beta, int tn)
// {
// 	double randgam_a,randgam_b;

// 	//error statment if either beta or alpha are <=0, as gamma distribution is undefined in this case
// 	if((beta<=0)||(alpha<=0))
// 	{
// 		ERR_CRITICAL("Beta distribution parameters in undefined range!\n");
// 	}

// 	randgam_a=gen_gamma_mt(1,alpha,tn);
// 	randgam_b=gen_gamma_mt(1,beta,tn);

// 	return randgam_a/(randgam_a+randgam_b);
// }

// /* function gen_lognormal(double mu, double sigma)
//  * purpose: to generate samples from a lognormal distribution with parameters mu and sigma
//  *
//  * parameters:
//  *  mean mu
//  *  standard deviation sigma
//  *
//  * returns: double from the specified lognormal distribution
//  *
//  * author: ggilani, date: 09/02/17
//  */
// double gen_lognormal(double mu,double sigma)
// {
// 	double randnorm,location,scale;

// 	randnorm=snorm();
// 	location=log(mu/sqrt(1+((sigma*sigma)/(mu*mu))));
// 	scale=sqrt(log(1+((sigma*sigma)/(mu*mu))));
// 	return exp(location+scale*randnorm);
// }

// /* function gen_lognormal(double mu, double sigma, int tn)
//  * purpose: to generate samples from a lognormal distribution with parameters mu and sigma
//  *
//  * parameters:
//  *  mean mu
//  *  standard deviation sigma
//  *  thread number tn
//  *
//  * returns: double from the specified lognormal distribution
//  *
//  * author: ggilani, date: 09/02/17
//  */
// double gen_lognormal_mt(double mu,double sigma,int tn)
// {
// 	double randnorm;

// 	randnorm=snorm_mt(tn);
// 	return exp(mu+sigma*randnorm);
// }

// void SampleWithoutReplacement(int tn, int k, int n)
// {
// /* Based on algorithm SG of http://portal.acm.org/citation.cfm?id=214402 
// ACM Transactions on Mathematical Software (TOMS) archive
// Volume 11 ,  Issue 2  (June 1985) table of contents
// Pages: 157 - 169   
// Year of Publication: 1985 
// ISSN:0098-3500 
// */

// 	double t,r,a,mu,f,es;
// 	int i,j,q,b,i2;

// 	if(k<3)
// 		{
// 		for(i=0;i<k;i++)
// 			{
// 			do
// 				{
// 				SamplingQueue[tn][i]=(int) (ranf_mt(tn)*((double) n));
// 				for(j=q=0;(j<i)&&(!q);j++)
// 					q=(SamplingQueue[tn][i]==SamplingQueue[tn][j]);
// 				}
// 			while(q);
// 			}
// 		q=k;
// 		}
// 	else if(2*k>n)
// 		{
// 		for(i=0;i<n;i++)
// 			SamplingQueue[tn][i]=i;
// 		for(i=n;i>k;i--)
// 			{
// 			j=(int) (ranf_mt(tn)*((double) i));
// 			if(j!=i-1)
// 				{
// 				b=SamplingQueue[tn][j];
// 				SamplingQueue[tn][j]=SamplingQueue[tn][i-1];
// 				SamplingQueue[tn][i-1]=b;
// 				}
// 			}
// 		q=k;
// 		}
// 	else if(4*k>n)
// 		{
// 		for(i=0;i<n;i++)
// 			SamplingQueue[tn][i]=i;
// 		for(i=0;i<k;i++)
// 			{
// 			j=(int) (ranf_mt(tn)*((double) (n-i)));
// 			if(j>0)
// 				{
// 				b=SamplingQueue[tn][i];
// 				SamplingQueue[tn][i]=SamplingQueue[tn][i+j];
// 				SamplingQueue[tn][i+j]=b;
// 				}
// 			}
// 		q=k;
// 		}
// 	else
// 		{
// 		/* fprintf(stderr,"@%i %i:",k,n); */
// 		t=(double) k;
// 		r=sqrt(t);
// 		a=sqrt(log(1+t/2*PI));
// 		a=a+a*a/(3*r);
// 		mu=t+a*r;
// 		b=2*MAX_PLACE_SIZE; /* (int) (k+4*a*r); */
// 		f=-1/(log(1-mu/((double) n)));
// 		i=q=0;
// 		while(i<=n) 
// 			{
// 			i+=(int) ceil(-log(ranf_mt(tn))*f);
// 			if(i<=n)
// 				{
// 				SamplingQueue[tn][q]=i-1;
// 				q++;
// 				if(q>=b) i=q=0;
// 				}
// 			else if(q<k)
// 				i=q=0;
// 			}
// 		}
// /*	else
// 		{
// 		t=(double) (n-k);
// 		r=sqrt(t);
// 		a=sqrt(log(1+t/2*PI));
// 		a=a+a*a/(3*r);
// 		mu=t+a*r;
// 		b=2*MAX_PLACE_SIZE; 
// 		f=-1/(log(1-mu/((double) n)));
// 		i=q=0;
// 		while(i<=n) 
// 			{
// 			i2=i+(int) ceil(-log(ranf_mt(tn))*f);
// 			i++;
// 			if(i2<=n)
// 				for(;(i<i2)&&(q<b);i++)
// 					{
// 					SamplingQueue[tn][q]=i-1;
// 					q++;
// 					}
// 			else 
// 				{
// 				for(;(i<=n)&&(q<b);i++)
// 					{
// 					SamplingQueue[tn][q]=i-1;
// 					q++;
// 					}
// 				if(q<k) i=q=0;
// 				}
// 			if(q>=b) i=q=0;
// 			}
// 		}
// */
// /*	if(k>2)
// 		{
// 		fprintf(stderr,"(%i) ",q);
// 		for(i=0;i<q;i++) fprintf(stderr,"%i ",SamplingQueue[tn][i]);
// 		fprintf(stderr,"\n");
// 		}
// */	while(q>k)
// 		{
// 		i=(int) (ranf_mt(tn)*((double) q));
// 		if(i<q-1) SamplingQueue[tn][i]=SamplingQueue[tn][q-1];
// 		q--;
// 		}

// }


/* RANDLIB static variables */
// long *Xcg1,*Xcg2;
// 	if(!(Xcg1=(long *) malloc(MAX_NUM_THREADS*CACHE_LINE_SIZE*sizeof(long)))) ERR_CRITICAL("Unable to allocate ranf storage\n");
// 	if(!(Xcg2=(long *) malloc(MAX_NUM_THREADS*CACHE_LINE_SIZE*sizeof(long)))) ERR_CRITICAL("Unable to allocate ranf storage\n");
// 	setall(P.seed1,P.seed2);

//Added Weibull pdf for calculating prob
// double WeibullCDF(int x, double shape, double scale)
// {
//     return (x < 0) ? 0 : (1-exp(-pow(((double)x/shape),scale)));
// }
