// Main header file
// Author: Yi Wang (yi dot wang at unsw dot edu dot au
// 16-Jun-2011

#ifndef _RESAMPTEST_H
#define _RESAMPTEST_H

// include R headers, but under 'strict' definitions
//
// RcppGSL includes GSL vector and matrix headers, plus
// Rcpp headers and a number of C++ standard library headers
#define STRICT_R_HEADERS
#include <RcppGSL.h>

#define printf Rprintf

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>

#include <sys/time.h>
#include <time.h>

// rmv.h
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_gamma.h>

// return status
#define SUCCESS 0
#define FAILED 1
#define CannotOpenFile 2
// logic
#define TRUE 1
#define FALSE 0
// model
#define LM 0
#define POISSON 1
#define NBIN 2
#define BIN 3
#define GAMMA 4
// link function
#define LOGIT 0
#define CLOGLOG 1
// shrinkage
#define NOSHRINK 0
#define IDENTITY 1
#define SHRINK 2
// test
#define LOGWILK 0
#define HOTELING 1
#define WALD 2
#define SCORE 3
#define LR 4
// estiMethod
#define NEWTON 0
#define CHI2 1
#define PHI 2
#define MOMENTS 3
// infoMatrix
#define OIM 0
#define EIM 1
// resampling
#define CASEBOOT 0
#define RESIBOOT 1
#define SCOREBOOT 2
#define PERMUTE 3
#define FREEPERM 4
#define MONTECARLO 5
#define RESIZ 6
#define SCOREZ 7
#define PITSBOOT 8
// p-value adjustment
#define NONE 0
#define UNADJUST 1
#define FREESTEP 2
#define SINGLESTEP 3
#define STEPUP 4
#define NOMONO 5
// R-squared
#define HOOPER 0
#define VECTOR 1
// others
#define TOL 1e-6
#define MAXITER 999
#define LAMBDA 0.8 // no shrinkage
#define NaN -100000
#define MAX_LINE_LENGTH 65536
#define WRAP 4
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(a) GSL_MAX(a, -a)
// simulations
#define NSIMU 500
#define NQ 6
#define RHO 0.5
#define ALFA 0.05
#define MAXLINELEN 256

typedef struct MethodStruc {
  // hypo test methods
  unsigned int nboot;
  unsigned int corr;
  unsigned int test;
  unsigned int resamp;
  unsigned int reprand;
  unsigned int student;
  unsigned int punit;
  unsigned int rsquare;
  unsigned int nRows;
  unsigned int nVars;
  unsigned int nParam;
  unsigned int showtime;
  unsigned int warning;

  // numeric
  double shrink_param;
  gsl_vector *smry_lambda;
  gsl_vector *anova_lambda;
  double tol;
} mv_Method;

// used for manylm only
typedef struct matStruc {
  gsl_matrix *mat; // hat(X)
  gsl_matrix *SS;
  gsl_matrix *R;
  gsl_matrix *Coef;
  gsl_matrix *Res;
  gsl_matrix *X;
  gsl_matrix *Y;
  gsl_vector *sd;
  double teststat;
} mv_mat;

typedef struct GroupMatrix {
  gsl_matrix *matrix;
} GrpMat;

// ManyGlm related
typedef struct RegressionMethod {
  // regression methods
  unsigned int model;
  unsigned int speclink;
  unsigned int varStab;
  unsigned int estiMethod;
  double tol;
  unsigned int maxiter, maxiter2;
  unsigned int n; // used in binomial regression
  unsigned int warning;
} reg_Method;

// ManyLM related
// anova.cpp
class AnovaTest {
public:
  mv_Method *mmRef;
  gsl_matrix *Yref;
  gsl_matrix *Xref;
  gsl_matrix *inRef;
  unsigned int nSamp;

  double *multstat;
  double *Pmultstat;
  gsl_matrix *statj;
  gsl_matrix *Pstatj;
  unsigned int *dfDiff;
  gsl_matrix *bootID;

  gsl_rng *rnd;
  // Methods
  AnovaTest(mv_Method *, gsl_matrix *, gsl_matrix *, gsl_matrix *isXvarIn);
  virtual ~AnovaTest();
  int resampTest(void);
  void releaseTest(void);
  void display(void);

private:
  mv_mat *Hats;
  gsl_permutation **sortid;
  gsl_vector *bStatj;
  double bMultStat;
  unsigned int nModels, nRows, nVars, nParam;

  // Methods
  int anovacase(gsl_matrix *bY, gsl_matrix *bX);
  int anovaresi(gsl_matrix *bY, const unsigned int p);
};

// summary.cpp
class Summary {
public:
  mv_Method *mmRef;
  gsl_matrix *Yref;
  gsl_matrix *Xref;
  unsigned int nSamp;
  double R2;
  double *multstat;
  double *Pmultstat;
  gsl_matrix *unitstat;
  gsl_matrix *Punitstat;
  gsl_matrix *bootID;
  gsl_rng *rnd;

  // Methods
  Summary(mv_Method *, gsl_matrix *, gsl_matrix *);
  virtual ~Summary();
  int resampTest(void);
  void releaseSummary(void);
  void display(void);

private:
  mv_mat *Hats;
  gsl_permutation **sortid;
  unsigned int nRows, nVars, nParam;
  double *bMultStat;
  gsl_matrix *bUnitStat;

  // Methods
  int calcR2(void);
  int smrycase(gsl_matrix *bY, gsl_matrix *bX);
  int smryresi(gsl_matrix *bY);
};

// glm base class
class glm {
public:
  glm(const reg_Method *mm);
  virtual ~glm();
  void initialGlm(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B);
  void releaseGlm(void);
  virtual int regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O,
                         gsl_matrix *B) = 0;
  virtual int EstIRLS(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O,
                      gsl_matrix *B, double *a) = 0;
  int copyGlm(glm *src);
  void display(void);

  // input arguments
  const reg_Method *mmRef;
  gsl_matrix *Yref;
  gsl_matrix *Xref;
  gsl_matrix *Oref;
  // return properties
  gsl_matrix *Beta;
  gsl_matrix *varBeta; // variance of Beta Hat
  gsl_matrix *Mu;
  gsl_matrix *Eta;
  gsl_matrix *Res;
  gsl_matrix *Var;
  gsl_matrix *wHalf;
  gsl_matrix *sqrt1_Hii;
  gsl_matrix *PitRes;

  unsigned int n;        // used in binomial regression
  unsigned int speclink; // used in binomial regression
  unsigned int rdf;
  double *theta, *ll, *dev, *aic;
  unsigned int *iterconv;
  unsigned int maxiter, maxiter2;
  double eps, mintol, maxtol, maxth;
  unsigned int nRows, nVars, nParams;
  // private:
  // abstract
  virtual double link(double) const = 0;
  virtual double invLink(double) const = 0;
  virtual double LinkDash(double) const = 0;
  virtual double weifunc(double, double) const = 0;
  virtual double varfunc(double, double) const = 0;
  virtual double llfunc(double, double, double) const = 0;
  virtual double devfunc(double, double, double) const = 0;
  virtual double pdf(double, double, double) const = 0;
  virtual double cdf(double, double, double) const = 0;
  virtual double cdfinv(double, double, double) const = 0;
  virtual double genRandist(double, double) const = 0;
};

// poisson regression
// this serves as a super class for the other families.
class PoissonGlm : public glm {
public: // public functions
  PoissonGlm(const reg_Method *);
  virtual ~PoissonGlm();
  virtual int regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O,
                         gsl_matrix *B) {
    return EstIRLS(Y, X, O, B, NULL);
  }
  int EstIRLS(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *, double *);
  int betaEst(unsigned int id, unsigned int iter, double *tol, double a);
  double getDisper(unsigned int id, double th) const;
  int update(gsl_vector *bj, unsigned int id);
  int predict(gsl_vector_view bj, unsigned int id, double a);
  // used by gamma
  double thetaEst_moments(unsigned int id);
  double thetaEst_newtons(double k0, unsigned int id, unsigned int limit);

  //    private:
  // Log-link and property functions
  double link(double mui) const { return log(mui); }
  double invLink(double etai) const { return exp(etai); }
  double LinkDash(double mui) const { return 1 / MAX(mui, mintol); }
  double weifunc(double mui, double a) const { return mui; }
  double varfunc(double mui, double a) const { return mui; }
  double llfunc(double yi, double mui, double a) const {
    if (yi == 0)
      return 2 * (-mui);
    else
      return 2 * (yi * log(mui) - mui - gsl_sf_lngamma(yi + 1));
  }
  double devfunc(double yi, double mui, double a) const {
    return 2 * (yi * log(MAX(1, yi) / mui) - yi + mui);
  }
  double pdf(double yi, double mui, double a) const {
    return Rf_dpois(yi, mui, FALSE);
  }
  double cdf(double yi, double mui, double a) const {
    return Rf_ppois(yi, mui, TRUE, FALSE);
  }
  double cdfinv(double u, double mui, double a) const {
    return Rf_qpois(u, mui, TRUE, FALSE);
  }
  double genRandist(double mui, double a) const { return Rf_rpois(mui); }
};

// Gamma regression y1 ~ GAMMA(shape = n = alpha = k, rate = shape / mui)
// Using log link
// See issue #53
class GammaGlm : public PoissonGlm {
public:
  GammaGlm(const reg_Method *);
  virtual ~GammaGlm();
  // we will use the log link inherited from PoissonGlm base class
  // the functions below are using the inverse link function
  // double link(double mui) const { return 1 / MAX(mintol, mui); }
  // double linkinv(double etai) const { return 1 / MAX(mintol, etai); }
  // double LinkDash(double mui) const { return -1 / MAX(mui * mui, mintol); }
  // weifunc is defined by the link not the family

  // The variance as a function of the mean
  // EX = shape / rate = mui, VarX = shape/(rate^2)
  // Var (mui) =  mui^2/shape
  // A here is the estimates shape parameter
  double varfunc(double mui, double a) const { return (mui * mui) / a; }
  // deviance residual HELP I am not sure if this is right.
  // see slide 28 of http://www.imm.dtu.dk/~hmad/GLM/slides/lect06.pdf
  double devfunc(double yi, double mui, double a) const {
    return -2 * (log(yi == 0 ? 1 : yi / MAX(mintol, mui)) -
                 (yi - mui) / MAX(mintol, mui));
  }
  // a = k = shape
  double llfunc(double yi, double mui, double a) const {
    return 2 * ((a - 1) * log(yi) - (yi * a) / mui +
                a * (log(a) - log(MAX(mintol, mui))) - gsl_sf_lngamma(a));
  }

  // mui = shape / rate => rate = shape / mui,
  // according to
  // http://dirk.eddelbuettel.com/code/rcpp/html/Rmath_8h_source.html Rf_ ...
  // takes the scale parameter and we use rate
  double pdf(double yi, double mui, double a) const {
    return Rf_dgamma(yi, a, mui / MAX(mintol, a), FALSE);
  }
  double cdf(double yi, double mui, double a) const {
    return Rf_pgamma(yi, a, mui / MAX(mintol, a), TRUE, FALSE);
  }
  double cdfinv(double ui, double mui, double a) const {
    return Rf_qgamma(ui, a, mui / MAX(mintol, a), TRUE, FALSE);
  }
  double genRandist(double mui, double a) const {
    return Rf_rgamma(a, mui / MAX(mintol, a));
  }
  // used by gamma
  double thetaEst_moments(unsigned int id);
  double thetaEst_newtons(double k0, unsigned int id, unsigned int limit);
};

// Binomial regression yi ~ BIN(n, pi), mui=n*pi
class BinGlm : public PoissonGlm {
public: // public functions
  BinGlm(const reg_Method *);
  virtual ~BinGlm();
  //    private: // property functions
  double link(double mui) const // pi=mui/n
  {
    if (speclink == CLOGLOG)
      return log(-log(1 - mui));
    else // default logit
      return log(mui / (n - mui));
  }
  double invLink(double ei) const {
    if (speclink == CLOGLOG)
      return MAX(mintol, MIN(1 - mintol, 1 - exp(-exp(ei))));
    else // default logit
      return n * exp(ei) / (1 + exp(ei));
  }
  double LinkDash(double mui) const {
    if (speclink == CLOGLOG)
      return 1 / MAX(mintol, (mui - 1) * log(1 - mui));
    else // default logit
      return n / MAX(mintol, mui * (n - mui));
  }
  double weifunc(double mui, double a) const {
    if (speclink == CLOGLOG)
      return ((1 - mui) * log(1 - mui)) * ((1 - mui) * log(1 - mui)) /
             MAX(mintol, (mui * (1 - mui)));
    else // default logit
      return mui * (1 - mui / n);
  }
  // others the same
  double varfunc(double mui, double a) const {
    return mui * (1 - mui / n);
  } // n*pi*(1-pi)
  double llfunc(double yi, double mui, double a) const {
    return 2 * ((yi > 0)
                    ? yi * log(mui / n)
                    : 0 + (yi < n) ? (n - yi) * log(1 - mui / n)
                                   : 0 + (n > 1) ? (gsl_sf_lngamma(n + 1) -
                                                    gsl_sf_lngamma(yi + 1) -
                                                    gsl_sf_lngamma(n - yi + 1))
                                                 : 0);
  }
  double devfunc(double yi, double mui, double a) const {
    return 2 * ((yi > 0) ? (yi * log(yi / mui))
                         : 0 + (yi < n) ? ((n - yi) * log((n - yi) / (n - mui)))
                                        : 0);
  }
  double pdf(double yi, double mui, double a) const
  //	        { if (n==1) return (yi<1)?(1-mui):mui;
  {
    return Rf_dbinom(yi, n, mui / n, FALSE);
  }
  double cdf(double yi, double mui, double a) const
  //                { if (n==1) return (yi<1)?(1-mui):1;
  {
    return Rf_pbinom(yi, n, mui / n, TRUE, FALSE);
  }
  double cdfinv(double ui, double mui, double a) const
  //                { if (n==1) return (ui<1)?0:1;
  {
    return (double)Rf_qbinom(ui, n, mui / n, TRUE, FALSE);
  }
  double genRandist(double mui, double a) const {
    return Rf_rbinom(n, mui / n);
  }
};

// negative binomial regression
class NBinGlm : public PoissonGlm // Y~NB(n, p)
{
public:
  NBinGlm(const reg_Method *mm);
  virtual ~NBinGlm();
  virtual int regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O,
                         gsl_matrix *B) {
    return nbinfit(Y, X, O, B);
  }
  int nbinfit(gsl_matrix *, gsl_matrix *, gsl_matrix *, gsl_matrix *);
  //    private:
  // See help(rnbinom) for the parameterisation used in ecology
  // *nbinom( , size, prob, mu, )
  // size = th (dispersion)
  // prob = size/(size+mu)
  double weifunc(double mui, double th) const {
    if (th == 0)
      return 0;
    else if (th > maxth)
      return mui; // poisson
    else
      return MAX(mintol, mui * th / MAX(mui + th, mintol));
  }
  double varfunc(double mui, double th) const {
    if (th == 0)
      return 0;
    else if (th > maxth)
      return mui;
    else
      return mui + mui * mui / MAX(th, mintol);
  }
  double llfunc(double yi, double mui, double th) const;
  double devfunc(double yi, double mui, double th) const {
    if (th == 0)
      return 0;
    else if (th > maxth)
      return 2 * (yi * log(MAX(1, yi) / mui) - yi + mui);
    else
      return 2 * (yi * log(MAX(1, yi) / mui) -
                  (yi + th) * log((yi + th) / (mui + th)));
  }
  double pdf(double yi, double mui, double th) const {
    if (th == 0)
      return 0;
    else if (th > maxth)
      return Rf_dpois(yi, mui, FALSE);
    else
      return Rf_dnbinom(yi, th, th / (mui + th), FALSE);
  }
  double cdf(double yi, double mui, double th) const {
    if (th == 0)
      return 1;
    else if (th > maxth)
      return Rf_ppois(yi, mui, TRUE, FALSE);
    else
      return Rf_pnbinom(yi, th, th / (mui + th), TRUE, FALSE);
  }
  double cdfinv(double ui, double mui, double th) const {
    if (th == 0)
      return 0;
    else if (th > maxth)
      return Rf_qpois(ui, mui, TRUE, FALSE);
    else
      return (double)Rf_qnbinom(ui, th, th / (mui + th), TRUE, FALSE);
  }
  double genRandist(double mui, double th) const {
    if (th == 0)
      return 0;
    else if (th > maxth)
      return Rf_rpois(mui);
    // Gamma-Poisson mixture, see rngbin()
    // rpois(lambda) where lambda=mu * rgamma(k, theta))/theta
    // for NB with mean mu and var=mu+mu^2/theta
    else
      return Rf_rnbinom(th, th / (mui + th));
  }

  double getfAfAdash(double th, unsigned int id, unsigned int limit);
  double thetaML(double th0, unsigned int id, unsigned int limit);
};

// base test class
class GlmTest {
public:
  const mv_Method *tm;
  glm *fit;
  gsl_matrix *Xin;
  gsl_matrix *smryStat, *Psmry;
  gsl_matrix *anovaStat, *Panova;
  gsl_matrix *bootID;

  // store the bootstrap stats
  gsl_matrix *bootStore;
  unsigned int nSamp;
  double *aic;
  unsigned int *dfDiff;

  // methods
  GlmTest(const mv_Method *tm);
  virtual ~GlmTest(void);
  void releaseTest(void);

  int summary(glm *);
  void displaySmry(glm *);

  int anova(glm *, gsl_matrix *);
  void displayAnova(void);

private:
  int getBootID(void);

  // int geeCalc(glm *PtrAlt, glm *PtrNull, gsl_matrix *);
  int GeeWald(glm *, gsl_matrix *, gsl_vector *);
  int GeeScore(gsl_matrix *, glm *, gsl_vector *);
  int GeeLR(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat);

  // int resampSmryCase(glm *, gsl_matrix *, GrpMat *, GrpMat *,
  // unsigned int i ); // summary
  int resampSmryCase(glm *, gsl_matrix *, GrpMat *, gsl_matrix *,
                     unsigned int i); // summary
  int resampAnovaCase(glm *, gsl_matrix *, gsl_matrix *, gsl_matrix *,
                      unsigned int i);
  int resampNonCase(glm *, gsl_matrix *, unsigned int i);

  // the following used in resampling
  unsigned int nModels;
  gsl_rng *rnd;
  size_t *permid;     // only useful in permutation test
  double lambda, eps; // intermediate shrinkage parameter

  // the following are used in geeCalc
  gsl_matrix *L; // only useful in Wald test
  gsl_matrix *Rlambda;
  gsl_matrix *Wj;
  gsl_matrix *XBeta, *Sigma; // used in monte carlo simulation
  // double *mr, *sr; // mean and variance of model residuals
  //	    gsl_vector *mr;

  GrpMat *GrpXs; // group X0
};

// io.cpp - input/output functions
int vector_filesize(FILE *f);
void matrix_filesize(FILE *f, int *row, int *col);
gsl_matrix *load_m(const char *file);
gsl_vector *load_v(const char *file);
void displaymatrix(gsl_matrix *m, const char *name);
void displayvector(gsl_vector *v, const char *name);

// calctest.c - utility functions
double logDet(gsl_matrix *SS);
int is_sym_matrix(const gsl_matrix *mat);
int subX(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subX1(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subX2(gsl_matrix *X, unsigned int id, gsl_matrix *Xi);
int subXrow(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int addXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subXrow1(gsl_matrix *X, gsl_vector *ref0, gsl_vector *ref1, gsl_matrix *Xi);
int GetR(gsl_matrix *Res, unsigned int corr, double lambda, gsl_matrix *R);
int subtractMean(gsl_matrix *dat);
// calctest.c - manylm related functions
int testStatCalc(mv_mat *H0, mv_mat *H1, mv_Method *mmRef,
                 const unsigned int ifcalcH1det, double *stat,
                 gsl_vector *statj);
int calcSS(gsl_matrix *Y, mv_mat *Hat, mv_Method *mmRef);
int calcAdjustP(const unsigned int punit, const unsigned int nVars, double *bj,
                double *sj, double *pj, gsl_permutation *sortid);
int reinforceP(double *p, unsigned int nVars, gsl_permutation *sortid);
int getHat(gsl_matrix *X, gsl_matrix *W, gsl_matrix *Hat);
int invLSQ(gsl_matrix *A, gsl_vector *b, gsl_vector *x);
int rcalc(gsl_matrix *Res, double, unsigned int, gsl_matrix *SS);

// simutility.cpp - functions used in simulation tests
int GetMean(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *Mu);
int GetPdstbtion(double *p, unsigned int nVars, unsigned int *isH0var,
                 unsigned int *cnt, unsigned int *cntfwe);
// int GetCov (gsl_matrix *Mu, gsl_matrix *Y, unsigned int AR1MAT, gsl_matrix
// *Sigma);  int GetMeanCov(gsl_matrix *X, gsl_matrix *Y, mv_Method *mm,
// unsigned int AR1MAT, gsl_matrix *Mu, gsl_matrix *Sigma);
int setMonteCarlo(glm *model, gsl_matrix *XBeta, gsl_matrix *Sigma);
int McSample(glm *model, gsl_rng *rnd, gsl_matrix *XBeta, gsl_matrix *Sigma,
             gsl_matrix *bT);

// rnd.c - functions to generate random numbers from multivariate (normal)
// distributions MVN random number generator
int rmvnorm(const gsl_rng *, const unsigned int, const gsl_matrix *,
            gsl_vector *);
// MVN with positive-semi definite covariance matrix
int semirmvnorm(const gsl_rng *, const unsigned int, const gsl_matrix *,
                gsl_vector *);

#endif
