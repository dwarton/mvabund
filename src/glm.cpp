// GLM estimation
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 16-333-2011

#include "resampTest.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <string.h>

// Note try to use gsl functions as much as possible to increase speed and
// stabilitiy

glm::glm(const reg_Method *mm)
    : mmRef(mm), Yref(NULL), Xref(NULL), Oref(NULL), Beta(NULL), varBeta(NULL),
      Mu(NULL), Eta(NULL), Res(NULL), Var(NULL), wHalf(NULL), sqrt1_Hii(NULL),
      PitRes(NULL), theta(NULL), ll(NULL), dev(NULL), aic(NULL),
      iterconv(NULL) {
  maxth = 100;
  n = mmRef->n; // the number of trials for a binomial repsonse .....
  speclink = mmRef->speclink;
  maxiter = mmRef->maxiter;
  maxiter2 = mmRef->maxiter2;
  // Error terms
  eps = mmRef->tol;
  // Valid data range of mu or y
  // mintol = mmRef->tol;
  mintol = 1e-6;
  if (mmRef->model == BIN)
    maxtol = n - mintol;
  else
    maxtol = 1 / mintol;
}

PoissonGlm::PoissonGlm(const reg_Method *mm) : glm(mm) {}

BinGlm::BinGlm(const reg_Method *mm) : PoissonGlm(mm) {}
GammaGlm::GammaGlm(const reg_Method *mm) : PoissonGlm(mm) {}

NBinGlm::NBinGlm(const reg_Method *mm) : PoissonGlm(mm) {}

glm::~glm() {}

PoissonGlm::~PoissonGlm() {}
GammaGlm::~GammaGlm() {}
BinGlm::~BinGlm() {}
NBinGlm::~NBinGlm() {}

void glm::releaseGlm(void) {
  if (Xref != NULL)
    gsl_matrix_free(Xref);
  if (Yref != NULL)
    gsl_matrix_free(Yref);
  if (Oref != NULL)
    gsl_matrix_free(Oref);
  if (Beta != NULL)
    gsl_matrix_free(Beta);
  if (varBeta != NULL)
    gsl_matrix_free(varBeta);
  if (Mu != NULL)
    gsl_matrix_free(Mu);
  if (Eta != NULL)
    gsl_matrix_free(Eta);
  if (Res != NULL)
    gsl_matrix_free(Res);
  if (Var != NULL)
    gsl_matrix_free(Var);
  if (wHalf != NULL)
    gsl_matrix_free(wHalf);
  if (sqrt1_Hii != NULL)
    gsl_matrix_free(sqrt1_Hii);
  if (PitRes != NULL)
    gsl_matrix_free(PitRes);
  if (theta != NULL)
    delete[] theta;
  if (ll != NULL)
    delete[] ll;
  if (dev != NULL)
    delete[] dev;
  if (iterconv != NULL)
    delete[] iterconv;
  if (aic != NULL)
    delete[] aic;
}

void glm::initialGlm(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O,
                     gsl_matrix *B) {
  releaseGlm();

  nRows = Y->size1;
  nVars = Y->size2;
  nParams = X->size2;

  unsigned int i, j;
  theta = new double[nVars];
  ll = new double[nVars];
  dev = new double[nVars];
  aic = new double[nVars];
  iterconv = new unsigned int[nVars];

  Xref = gsl_matrix_alloc(nRows, nParams);
  gsl_matrix_memcpy(Xref, X);
  Yref = gsl_matrix_alloc(nRows, nVars);
  gsl_matrix_memcpy(Yref, Y);
  if (O == NULL)
    Oref = NULL;
  else {
    Oref = gsl_matrix_alloc(nRows, nVars);
    gsl_matrix_memcpy(Oref, O);
  }
  Beta = gsl_matrix_alloc(nParams, nVars);
  varBeta = gsl_matrix_alloc(nParams, nVars);
  Mu = gsl_matrix_alloc(nRows, nVars);
  Eta = gsl_matrix_alloc(nRows, nVars);
  Res = gsl_matrix_alloc(nRows, nVars);
  Var = gsl_matrix_alloc(nRows, nVars);
  wHalf = gsl_matrix_alloc(nRows, nVars);
  sqrt1_Hii = gsl_matrix_alloc(nRows, nVars);
  PitRes = gsl_matrix_alloc(nRows, nVars);

  gsl_matrix_set_zero(varBeta);

  for (j = 0; j < nVars; j++) {
    theta[j] = maxtol; // i.e. phi=0
    ll[j] = 0;
    dev[j] = 0;
    aic[j] = 0;
    iterconv[j] = 0;
  }
  //  Note: setting the initial value is important
  //  e.g., using mean(Y) for binomial regression doesn't work
  //    gsl_matrix *t1;
  //    t1 = gsl_matrix_alloc(nRows, 1);
  //    gsl_matrix_set_all (t1, 1.0); // intercept
  //    GetMean(t1, Y, Mu);
  //    gsl_matrix_free(t1);
  //
  //  Use binomial$initialize: MuStart = (Y+0.5)/2
  //  It seems to work for poisson and negative.binomial as well
  double LinAdjust, ScaleAdjust;
  double eij;
  // cLogLog initialization is the same as binomial
  if (mmRef->model == BIN) {
    LinAdjust = 0.5;
    ScaleAdjust = 0.5;
  } else {
    LinAdjust = 0.1;
    ScaleAdjust = 1;
  }
  if (B != NULL) {
    gsl_matrix_memcpy(Beta, B);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, X, Beta, 0.0, Eta);
    for (i = 0; i < nRows; i++)
      for (j = 0; j < nVars; j++) {
        eij = gsl_matrix_get(Eta, i, j);
        // to avoid nan
        if (eij > link(maxtol))
          eij = link(maxtol);
        if (eij < link(mintol))
          eij = link(mintol);
        gsl_matrix_set(Eta, i, j, eij);
        gsl_matrix_set(Mu, i, j, invLink(eij));
      }
  } else if (O != NULL) {
    gsl_matrix_set_zero(Beta);
    gsl_matrix_memcpy(Eta, O);
    for (i = 0; i < nRows; i++)
      for (j = 0; j < nVars; j++) {
        eij = gsl_matrix_get(Eta, i, j);
        // to avoid nan
        if (eij > link(maxtol))
          eij = link(maxtol);
        if (eij < link(mintol))
          eij = link(mintol);
        gsl_matrix_set(Eta, i, j, eij);
        gsl_matrix_set(Mu, i, j, invLink(eij));
      }
  } else {
    gsl_matrix_memcpy(Mu, Yref);
    gsl_matrix_add_constant(Mu, LinAdjust);
    gsl_matrix_scale(Mu, ScaleAdjust);
    //       gsl_matrix_set_zero (Eta); // intercept
    for (i = 0; i < nRows; i++)
      for (j = 0; j < nVars; j++) {
        eij = link(gsl_matrix_get(Mu, i, j));
        if (eij > link(maxtol))
          eij = link(maxtol);
        if (eij < link(mintol))
          eij = link(mintol);
        gsl_matrix_set(Eta, i, j, eij);
        gsl_matrix_set(Mu, i, j, invLink(eij));
      }
    gsl_matrix_set_zero(Beta); // intercept
    gsl_vector_view b0 = gsl_matrix_column(Beta, 0);
    gsl_vector_set_all(&b0.vector, 1.0);
    //    printf("LinAdjust=%.2f, ScaleAdjust=%.2f\n", LinAdjust, ScaleAdjust);
  }

  rdf = nRows - nParams;
}

int glm::copyGlm(glm *src) {
  initialGlm(src->Yref, src->Xref, src->Oref, NULL);

  // copy properties
  gsl_matrix_memcpy(Beta, src->Beta);
  gsl_matrix_memcpy(varBeta, src->varBeta);
  gsl_matrix_memcpy(Mu, src->Mu);
  gsl_matrix_memcpy(Eta, src->Eta);
  gsl_matrix_memcpy(Res, src->Res);
  gsl_matrix_memcpy(Var, src->Var);
  gsl_matrix_memcpy(wHalf, src->wHalf);
  gsl_matrix_memcpy(sqrt1_Hii, src->sqrt1_Hii);
  gsl_matrix_memcpy(PitRes, src->PitRes);

  for (unsigned int i = 0; i < nVars; i++) {
    theta[i] = src->theta[i];
    ll[i] = src->ll[i];
    dev[i] = src->dev[i];
    iterconv[i] = src->iterconv[i];
    aic[i] = src->aic[i];
  }

  return SUCCESS;
}

int PoissonGlm::EstIRLS(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O,
                        gsl_matrix *B, double *a) {
  initialGlm(Y, X, O, B);

  gsl_set_error_handler_off();
  gsl_rng *rnd = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int i, j;
  int status;
  double yij, mij, vij, wij, tol, hii, uij, wei;
  gsl_vector_view Xwi, Xi, vj, hj, dj;

  gsl_matrix *WX = gsl_matrix_alloc(nRows, nParams);
  gsl_matrix *TMP = gsl_matrix_alloc(nRows, nParams);
  gsl_matrix *XwX = gsl_matrix_alloc(nParams, nParams);
  for (j = 0; j < nVars; j++) {
    if (a != NULL)
      theta[j] = a[j];
    // estimate mu and beta
    iterconv[j] = betaEst(j, maxiter, &tol, theta[j]);
    // estimate the shape for gamma
    if (mmRef->model == GAMMA) {
      theta[j] = mmRef->estiMethod == MOMENTS
                     ? thetaEst_moments(j)
                     : thetaEst_newtons(0, j, maxiter2);
    }
    if ((mmRef->warning == TRUE) & (iterconv[j] == maxiter))
      printf("Warning: EstIRLS reached max iterations, may not converge in the "
             "%d-th variable (dev=%.4f, err=%.4f)!\n",
             j, dev[j], tol);
    gsl_matrix_memcpy(WX, X);
    for (i = 0; i < nRows; i++) {
      mij = gsl_matrix_get(Mu, i, j);
      // get variance
      vij = varfunc(mij, theta[j]);
      gsl_matrix_set(Var, i, j, vij);
      // get weight
      wij = sqrt(weifunc(mij, theta[j]));
      gsl_matrix_set(wHalf, i, j, wij);
      // get (Pearson) residuals
      yij = gsl_matrix_get(Y, i, j);
      gsl_matrix_set(Res, i, j, (yij - mij) / sqrt(vij));
      // PIT residuals
      // get PIT residuals for cts families, ie just the cdf of the observations
      if (mmRef->model == GAMMA) {
        uij = cdf(yij, mij, theta[j]);
      } else {
        // get PIT residuals for discrete data
        wei = gsl_rng_uniform_pos(rnd); // wei ~ U(0, 1)
        uij = wei * cdf(yij, mij, theta[j]);
        if (yij > 0)
          uij = uij + (1 - wei) * cdf((yij - 1), mij, theta[j]);
      }
      gsl_matrix_set(PitRes, i, j, uij);
      // get elementry log-likelihood
      ll[j] = ll[j] + llfunc(yij, mij, theta[j]);
      // W^1/2 X
      Xwi = gsl_matrix_row(WX, i);
      gsl_vector_scale(&Xwi.vector, wij);
    }
    aic[j] = -ll[j] + 2 * (nParams);

    // X^T * W * X
    gsl_matrix_set_identity(XwX);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);
    status = gsl_linalg_cholesky_decomp(XwX);
    if (status == GSL_EDOM) {
      if (mmRef->warning == TRUE)
        printf("Warning: singular matrix in calculating pit-residuals. An "
               "eps*I is added to the singular matrix.\n");
      gsl_matrix_set_identity(XwX);
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, WX, mintol, XwX);
      gsl_linalg_cholesky_decomp(XwX);
    }
    gsl_linalg_cholesky_invert(XwX);

    // Calc varBeta
    dj = gsl_matrix_diagonal(XwX);
    vj = gsl_matrix_column(varBeta, j);
    gsl_vector_memcpy(&vj.vector, &dj.vector);

    // hii is diagonal element of H=X*(X'WX)^-1*X'*W
    hj = gsl_matrix_column(sqrt1_Hii, j);
    gsl_blas_dsymm(CblasRight, CblasLower, 1.0, XwX, Xref, 0.0,
                   TMP); // X*(X'WX)^-1
    for (i = 0; i < nRows; i++) {
      Xwi = gsl_matrix_row(TMP, i);
      Xi = gsl_matrix_row(Xref, i);
      wij = gsl_matrix_get(wHalf, i, j);
      gsl_blas_ddot(&Xwi.vector, &Xi.vector, &hii);
      gsl_vector_set(&hj.vector, i,
                     MAX(mintol, sqrt(MAX(0, 1 - wij * wij * hii))));
    }
  }
  // standardize perason residuals by rp/sqrt(1-hii)
  //   gsl_matrix_div_elements (Res, sqrt1_Hii);
  //   subtractMean(Res);  // have mean subtracted

  gsl_matrix_free(XwX);
  gsl_matrix_free(WX);
  gsl_matrix_free(TMP);
  gsl_rng_free(rnd);

  return SUCCESS;
}

int PoissonGlm::betaEst(unsigned int id, unsigned int iter, double *tol,
                        double th) {
  gsl_set_error_handler_off();
  int status, isValid;
  // unsigned int j, ngoodobs;
  unsigned int i, step, step1;
  double wij, zij, eij, mij, yij; //, bij;
  double dev_old, dev_grad = 1.0;
  gsl_vector_view Xwi;
  gsl_matrix *WX, *XwX;
  gsl_vector *z, *Xwz;
  gsl_vector *coef_old = gsl_vector_alloc(nParams);
  gsl_vector_view bj = gsl_matrix_column(Beta, id);

  // Main Loop of IRLS begins
  z = gsl_vector_alloc(nRows);
  WX = gsl_matrix_alloc(nRows, nParams);
  XwX = gsl_matrix_alloc(nParams, nParams);
  Xwz = gsl_vector_alloc(nParams);
  step = 0;
  *tol = 1.0;
  gsl_vector_memcpy(coef_old, &bj.vector);
  while (step < iter) {
    for (i = 0; i < nRows; i++) {
      // (y-m)/g'
      yij = gsl_matrix_get(Yref, i, id);
      eij = gsl_matrix_get(Eta, i, id);
      mij = gsl_matrix_get(Mu, i, id);
      // if (mij<mintol) mij=mintol;
      // if (mij>maxtol) mij=maxtol;
      zij = eij + (yij - mij) * LinkDash(mij);
      if (Oref != NULL) {
        zij = zij - gsl_matrix_get(Oref, i, id);
      }
      // wt=sqrt(weifunc);
      wij = sqrt(weifunc(mij, th));
      // W^1/2*z[good]
      gsl_vector_set(z, i, wij * zij);
      // W^1/2*X[good]
      Xwi = gsl_matrix_row(Xref, i);
      gsl_matrix_set_row(WX, i, &Xwi.vector);
      Xwi = gsl_matrix_row(WX, i);
      gsl_vector_scale(&Xwi.vector, wij);
    }
    // in glm2, solve WXb=Wz, David suggested not good
    // So back to solving X'WXb=X'Wz
    gsl_matrix_set_identity(XwX);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);
    status = gsl_linalg_cholesky_decomp(XwX);
    if (status == GSL_EDOM) {
      if (mmRef->warning == TRUE) {
        printf("Warning: singular matrix in betaEst: ");
        gsl_matrix_set_identity(XwX);
        gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, Xref, 0.0, XwX);
        //  displaymatrix(Xref, "Xref");
        // displaymatrix(XwX, "XX^T");
        //   printf("calc(XX')=%.8f\n", calcDet(XwX));
        status = gsl_linalg_cholesky_decomp(XwX);
        if (status == GSL_EDOM)
          printf("X^TX is singular - check case resampling or input design "
                 "matrix!\n");
        else {
          for (i = 0; i < nRows; i++) {
            mij = gsl_matrix_get(Mu, i, id);
            wij = sqrt(weifunc(mij, th));
            if (wij < mintol)
              printf("weight[%d, %d]=%.4f is too close to zero\n", i, id, wij);
          }
        }
        printf("An eps*I is added to the singular matrix.\n");
      }
      gsl_matrix_set_identity(XwX);
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, WX, mintol, XwX);
      gsl_linalg_cholesky_decomp(XwX);
    }
    gsl_blas_dgemv(CblasTrans, 1.0, WX, z, 0.0, Xwz);
    gsl_linalg_cholesky_solve(XwX, Xwz, &bj.vector);

    // Debug for nan
    /*
    if (gsl_vector_get(&bj.vector, 1) != gsl_vector_get(&bj.vector, 1)) {
      displayvector(z, "z");
      displayvector(&bj.vector, "bj");
      gsl_vector_view mj = gsl_matrix_column(Mu, id);
      displayvector(&mj.vector, "mj");
      printf("weight\n");
      for (i = 0; i < nRows; i++) {
        printf("%.4f ", sqrt(weifunc(mij, th)));
      }
      printf("\n");
      displaymatrix(XwX, "XwX");
      exit(-1);
    }
    */
    // Given bj, update eta, mu
    dev_old = dev[id];
    isValid = predict(bj, id, th);
    dev_grad = (dev[id] - dev_old) / (ABS(dev[id]) + 0.1);
    *(tol) = ABS(dev_grad);

    step1 = 0;
    // If divergent or increasing deviance, half step
    // (step>1) -> (step>0) gives weired results for NBin fit
    // below works for boundary values, esp BIN fit but not NBin fit
    while ((dev_grad > eps) & (step > 1)) {
      gsl_vector_add(&bj.vector, coef_old);
      gsl_vector_scale(&bj.vector, 0.5);
      //     dev_old=dev[id];
      isValid = predict(bj, id, th);
      dev_grad = (dev[id] - dev_old) / (ABS(dev[id]) + 0.1);
      *tol = ABS(dev_grad);
      if (*tol < eps)
        break;
      step1++;
      if (step1 > 10) {
        // printf("\t Half step stopped at iter %d: gradient=%.8f\n", step1,
        //       dev_grad);
        break;
      }
    }
    if (isValid == TRUE)
      gsl_vector_memcpy(coef_old, &bj.vector);

    step++;
    if (*tol < eps)
      break;
  }

  gsl_vector_free(z);
  gsl_matrix_free(WX);
  gsl_matrix_free(XwX);
  gsl_vector_free(Xwz);
  gsl_vector_free(coef_old);

  return step;
}

int PoissonGlm::update(gsl_vector *bj, unsigned int id) {
  int isValid = TRUE;
  unsigned int i;
  double eij, mij;
  gsl_vector_view xi;

  for (i = 0; i < nRows; i++) {
    xi = gsl_matrix_row(Xref, i);
    gsl_blas_ddot(&xi.vector, bj, &eij); // set eij to xi \dot bi
    if (Oref != NULL)
      eij = eij + gsl_matrix_get(Oref, i, id);
    if (eij > link(maxtol)) { // to avoid nan;
      eij = link(maxtol);
      isValid = FALSE;
    }
    if (eij < link(mintol)) {
      eij = link(mintol);
      isValid = FALSE;
    }
    mij = invLink(eij);
    gsl_matrix_set(Eta, i, id, eij);
    gsl_matrix_set(Mu, i, id, mij);
  }

  return isValid;
}

int PoissonGlm::predict(gsl_vector_view bj, unsigned int id, double th) {
  unsigned int i;
  double yij, mij;

  int isValid = update(&bj.vector, id);

  // Given valid mj, estimate deviance
  dev[id] = 0;
  for (i = 0; i < nRows; i++) {
    yij = gsl_matrix_get(Yref, i, id);
    mij = gsl_matrix_get(Mu, i, id);
    dev[id] = dev[id] + devfunc(yij, mij, th);

    // if (mmRef->model == GAMMA && mmRef->warning == TRUE &&
    //    i % 25 == 0) // REMOVE
    //  printf("dev = %f yij = %f, mij = %f\n", devfunc(yij, mij, th), yij, mij,
    //         th);
  }

  return isValid;
}

// estimates the shape (k)  parameter for the gamma family
// via the newtons method, starting with a method of moments estimate
// https://tminka.github.io/papers/minka-gamma.pdf
// Newtons method outlined here
// https://www.tandfonline.com/doi/abs/10.1080/00401706.1969.10490731
// Method of moments estimate from slide 22:
// https://www.statistics.ma.tum.de/fileadmin/w00bdb/www/czado/lec8.pdf
// k- is starting
double PoissonGlm::thetaEst_moments(unsigned int id) {
  unsigned int i;
  double sum = 0, num = 0;
  double y, m;

  for (i = 0; i < nRows; i++) {
    y = gsl_matrix_get(Yref, i, id);
    m = gsl_matrix_get(Mu, i, id);
    sum = sum + (y / m - 1) * (y / m - 1);
    num = num + 1;
  }
  return (num - nParams) / sum;
}
double PoissonGlm::thetaEst_newtons(double k0, unsigned int id,
                                    unsigned int limit) {
  unsigned int i, it = 0;
  double y;
  // obtain an initial estimate via MoM
  double k = k0 == 0 ? thetaEst_moments(id) : k0;
  // calculate some constants
  double logxbar = 0, barlogx = 0;
  for (i = 0; i < nRows; i++) {
    y = gsl_matrix_get(Yref, i, id);
    logxbar += y;
    barlogx += log(y);
  }
  logxbar = log(logxbar / nRows);
  barlogx = barlogx / nRows;
  double s = logxbar - barlogx;
  // now use newton-raphson
  while (it < limit) {
    it++;
    double update = (log(k) - gsl_sf_psi(k) - s) / ((1 / k) - gsl_sf_psi_1(k));
    k = k - update;
    double tol = ABS(update);
    // break if the update was very small
    if (tol < eps)
      break;
  }
  return k;
}

int NBinGlm::nbinfit(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O,
                     gsl_matrix *B) {
  gsl_set_error_handler_off();

  initialGlm(Y, X, O, B);

  gsl_rng *rnd = gsl_rng_alloc(gsl_rng_mt19937);
  unsigned int i, j; //, isConv;
  double yij, mij, vij, hii, uij, wij, wei;
  double th, tol, dev_th_b_old;
  int status;
  //   gsl_vector_view b0j, m0j, e0j, v0j;
  gsl_matrix *WX = gsl_matrix_alloc(nRows, nParams);
  gsl_matrix *TMP = gsl_matrix_alloc(nRows, nParams);
  gsl_matrix *XwX = gsl_matrix_alloc(nParams, nParams);
  gsl_vector_view Xwi, Xi, vj, dj, hj;

  for (j = 0; j < nVars; j++) {
    betaEst(j, maxiter, &tol, maxtol); // poisson
    // Get initial theta estimates
    iterconv[j] = 0.0;
    if (mmRef->estiMethod == CHI2) {
      th = 1.0; // initial phi value
      while (iterconv[j] < maxiter) {
        iterconv[j]++;
        dev_th_b_old = dev[j];
        th = th / getDisper(j, th); // Note: getDisper returns phi=1/theta!
        betaEst(j, 1.0, &tol, th);  // 1-step beta
        tol = ABS((dev[j] - dev_th_b_old) / (ABS(dev[j]) + 0.1));
        if (tol < eps)
          break;
      }
    } else if (mmRef->estiMethod == NEWTON) {
      th = thetaML(0.0, j, maxiter);
      while (iterconv[j] < maxiter) {
        iterconv[j]++;
        dev_th_b_old = dev[j];
        th = thetaML(th, j, maxiter2);
        betaEst(j, maxiter2, &tol, th);
        tol = ABS((dev[j] - dev_th_b_old) / (ABS(dev[j]) + 0.1));
        if (tol < eps)
          break;
      }
    } else {
      th = getfAfAdash(0.0, j, maxiter);
      /*           lm=0;
                 for (i=0; i<nRows; i++) {
                     yij = gsl_matrix_get(Y, i, j);
                     mij = gsl_matrix_get(Mu, i, j);
                     lm = lm + llfunc( yij, mij, th);
                 } */
      while (iterconv[j] < maxiter) {
        iterconv[j]++;
        dev_th_b_old = dev[j];
        betaEst(j, maxiter2, &tol, th);
        th = getfAfAdash(th, j, 1.0);
        tol = ABS((dev[j] - dev_th_b_old) / (ABS(dev[j]) + 0.1));
        if (tol < eps)
          break;
      }
    }
    if ((iterconv[j] == maxiter) & (mmRef->warning == TRUE))
      printf("Warning: reached maximum itrations - negative binomial may NOT "
             "converge in the %d-th variable (dev=%.4f, err=%.4f, "
             "theta=%.4f)!\n",
             j, dev[j], tol, th);

    // other properties based on mu and phi
    theta[j] = th;
    gsl_matrix_memcpy(WX, Xref);
    ll[j] = 0;
    for (i = 0; i < nRows; i++) {
      yij = gsl_matrix_get(Y, i, j);
      mij = gsl_matrix_get(Mu, i, j);
      vij = varfunc(mij, th);
      gsl_matrix_set(Var, i, j, vij);
      wij = sqrt(weifunc(mij, th));
      gsl_matrix_set(wHalf, i, j, wij);
      gsl_matrix_set(Res, i, j, (yij - mij) / sqrt(vij));
      ll[j] = ll[j] + llfunc(yij, mij, th);
      // get PIT residuals for discrete data
      wei = gsl_rng_uniform_pos(rnd); // wei ~ U(0, 1)
      uij = wei * cdf(yij, mij, th);
      if (yij > 0)
        uij = uij + (1 - wei) * cdf((yij - 1), mij, th);
      gsl_matrix_set(PitRes, i, j, uij);
      // W^1/2 X
      Xwi = gsl_matrix_row(WX, i);
      gsl_vector_scale(&Xwi.vector, wij);
    }
    aic[j] = -ll[j] + 2 * (nParams + 1);

    // X^T * W * X
    gsl_matrix_set_identity(XwX);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);
    status = gsl_linalg_cholesky_decomp(XwX);
    if (status == GSL_EDOM) {
      if (mmRef->warning == TRUE)
        printf("Warning: singular matrix in calculating pit-residuals. An "
               "eps*I is added to the singular matrix.\n");
      gsl_matrix_set_identity(XwX);
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, WX, mintol, XwX);
      gsl_linalg_cholesky_decomp(XwX);
    }
    gsl_linalg_cholesky_invert(XwX); // (X'WX)^-1

    // Calc varBeta
    vj = gsl_matrix_column(varBeta, j);
    dj = gsl_matrix_diagonal(XwX);
    gsl_vector_memcpy(&vj.vector, &dj.vector);

    // hii is diagonal element of H=X*(X'WX)^-1*X'*W
    hj = gsl_matrix_column(sqrt1_Hii, j);
    gsl_blas_dsymm(CblasRight, CblasLower, 1.0, XwX, Xref, 0.0,
                   TMP); // X*(X'WX)^-1
    for (i = 0; i < nRows; i++) {
      Xwi = gsl_matrix_row(TMP, i);
      Xi = gsl_matrix_row(Xref, i);
      wij = gsl_matrix_get(wHalf, i, j);
      gsl_blas_ddot(&Xwi.vector, &Xi.vector, &hii);
      gsl_vector_set(&hj.vector, i,
                     MAX(mintol, sqrt(MAX(0, 1 - wij * wij * hii))));
      // printf("hii=%.4f, wij=%.4f, sqrt(1-wij*wij*hii)=%.4f\n", hii, wij,
      // sqrt(1-wij*wij*hii));
    }
  } // end nVar for j loop
    //   gsl_matrix_div_elements (Res, sqrt1_Hii);
    //   subtractMean(Res);

  gsl_matrix_free(XwX);
  gsl_matrix_free(WX);
  gsl_matrix_free(TMP);
  gsl_rng_free(rnd);

  return SUCCESS;
}

// old based on the phi-1/theta on log scale in <Hilbe2008>
double PoissonGlm::getDisper(unsigned int id, double th) const {
  unsigned int i, df, nNonZero = 0;
  double ss2, yij, mij, chi2 = 0, var = 1.0;

  gsl_vector_view yj = gsl_matrix_column(Yref, id);
  gsl_vector_view mj = gsl_matrix_column(Mu, id);
  for (i = 0; i < nRows; i++) {
    yij = gsl_vector_get(&yj.vector, i);
    mij = gsl_vector_get(&mj.vector, i);
    ss2 = (yij - mij) * (yij - mij); // ss = (y-mu)^2
    if (mij < mintol)
      mij = mintol; // mij = 1;
    else
      nNonZero++;
    if (var > eps) {
      var = varfunc(mij, th);
      chi2 = chi2 + ss2 / var; // dist dependant
    }
  }
  if (nNonZero > nParams)
    df = nNonZero - nParams;
  else
    df = 1;
  //    df = nRows - nParams;
  return chi2 / df;
}

double NBinGlm::llfunc(double yi, double mui, double th) const {
  double l = 0; //, p=0;

  if (th == 0) {
    l = gsl_sf_lngamma(mintol) + log(MAX(yi, mintol));
  } else if (th > maxth) { // Poisson
    l = -yi * log(mui) + mui + gsl_sf_lngamma(yi + 1);
  } else {
    l = (yi + th) * log(mui + th) - yi * log(mui) + gsl_sf_lngamma(yi + 1) -
        th * log(th) + gsl_sf_lngamma(th) - gsl_sf_lngamma(yi + th);
  }

  if (l != l)
    printf("l=nan, theta=%.4f, yi=%.4f, mu=%.4f\n", th, yi, mui);

  return -2 * l;
}

double NBinGlm::thetaML(double k0, unsigned int id, unsigned int limit) {
  // equivalent to theta.ml() in MASS
  // Note that theta here is the dispersion parameter
  // So phi = 1/theta;
  unsigned int i, it = 0;
  double del = 1, sum = 1, num = 0, k;
  double y, m, dl, ddl, tol;
  // inital guess via method of moments (?maybe)
  if (k0 == 0) {
    for (i = 0; i < nRows; i++) {
      y = gsl_matrix_get(Yref, i, id);
      m = gsl_matrix_get(Mu, i, id);
      if (m > 0) {
        sum = sum + (y / m - 1) * (y / m - 1);
        num = num + 1;
      }
    }
    k = num / sum;
  } else
    k = k0;

  k = MAX(k, mintol);
  while (it <= limit) {
    it++;
    k = ABS(k);
    dl = nRows * (1 + log(k) - gsl_sf_psi(k));
    ddl = nRows * (gsl_sf_psi_1(k) - 1 / k);
    for (i = 0; i < nRows; i++) {
      y = gsl_matrix_get(Yref, i, id);
      m = gsl_matrix_get(Mu, i, id);
      dl = dl + gsl_sf_psi(y + k) - log(m + k) - (y + k) / (m + k);
      ddl = ddl - gsl_sf_psi_1(y + k) + 2 / (m + k) -
            (y + k) / ((m + k) * (m + k));
    }
    if (ABS(ddl) < mintol)
      ddl = GSL_SIGN(ddl) * mintol;
    del = dl / ABS(ddl);
    tol = ABS(del * dl);
    if (tol < eps)
      break;
    k = k + del; // Normal Newton use - instead of + for -ddl
    if (k > maxth)
      break;
    if (k < 0) {
      k = 0;
      break;
    }
  }
  // if (k<0) k=0;

  return k;
}

// moments estimation (maybe)
double NBinGlm::getfAfAdash(double k0, unsigned int id, unsigned int limit) {
  unsigned int i, it = 0;
  double sum = 1, num = 0, k;
  double y, m, dl, ddl, tol;
  double phi, dl_dphi, d2l_dphi2, del_phi;
  if (k0 == 0) {
    for (i = 0; i < nRows; i++) {
      y = gsl_matrix_get(Yref, i, id);
      m = gsl_matrix_get(Mu, i, id);
      if (m > 0) {
        sum = sum + (y / m - 1) * (y / m - 1);
        num = num + 1;
      }
    }
    k = num / sum;
    if (num == 0)
      printf("num=0\n");
  } else
    k = k0;
  k = MAX(k, mintol);
  phi = 1 / k;
  while (it < limit) {
    it++;
    dl = nRows * (1 + log(k) - gsl_sf_psi(k));
    ddl = nRows * (gsl_sf_psi_1(k) - 1 / k);
    for (i = 0; i < nRows; i++) {
      y = gsl_matrix_get(Yref, i, id);
      m = gsl_matrix_get(Mu, i, id);
      dl = dl + gsl_sf_psi(y + k) - log(m + k) - (y + k) / (m + k);
      ddl = ddl - gsl_sf_psi_1(y + k) + 2 / (m + k) -
            (y + k) / ((m + k) * (m + k));
    }
    dl_dphi = -exp(2 * log(k)) * dl;
    d2l_dphi2 = 2 * exp(3 * log(k)) * dl + exp(4 * log(k)) * ddl;

    if (ABS(ddl) < mintol)
      ddl = GSL_SIGN(ddl) * mintol;
    del_phi = dl_dphi / ABS(d2l_dphi2);
    tol = ABS(del_phi * dl_dphi);

    if (tol < eps)
      break;

    phi = phi + del_phi;
    if (phi < 0) {
      k = 0;
      break;
    }
    k = 1 / MAX(ABS(phi), mintol);
    if (k > maxth)
      break;
  }

  return k;
}

void glm::display(void) {
  unsigned int j;
  /*
      if ( mmRef->model == LM )
         printf("Linear regression:\n");
      else if ( mmRef->model == POISSON )
         printf("Poisson regression:\n");
      else if ( mmRef->model == BIN ) {
         if (n==1 ) printf("Logistic regression:\n");
         else printf("Binomial regression:\n");
      }
      else if ( mmRef->model == NB ) {
         printf("Negative Binomial regression ");
         switch (mmRef->estiMethod) {
            case NEWTON:
                printf("(Newton-ML):\n");
            break;
        case CHI2:
            printf("(Chi2):\n");
            break;
            case FISHER:
                printf("(Fisher Scoring):\n");
            break;
            default:
                printf("phi estimation method not available");
         }
      }
      else {
          printf("GLM regression method not available");
      }
  */
  printf("Two-log-like=\n ");
  for (j = 0; j < nVars; j++)
    printf("%.2f ", ll[j]);
  printf("\n");
  //    printf("AIC=\n " );
  //    for ( j=0; j<nVars; j++ ) printf("%.2f ", aic[j]);
  //    printf("\n");
  //    printf("# of convergence\n");
  //    for ( j=0; j<nVars; j++ )
  //        printf("%d ", iterconv[j]);
  //    printf("\n");
  //   printf("Residual deviance=\n " );
  //    for ( j=0; j<nVars; j++ ) printf("%.2f ", dev[j]);
  //    printf("\n");
  //    if ( mmRef->model == NB ) {
  //        printf("\nphi=\n ");
  //        for (j=0; j<nVars; j++ ) printf("%.2f ", phi[j]);
  //    }
  //    printf("\n");
  //    if (Oref != NULL) displaymatrix(Oref, "O");
  //    displaymatrix(Xref, "X");
  //    displaymatrix(Eta, "Eta");
  //    displaymatrix(Beta, "Beta");
  //    displaymatrix(varBeta, "varBeta");
  //    displaymatrix(Mu, "Mu");
  //    displaymatrix(Var, "Var");
  //    displaymatrix(PitRes, "PitRes");
  //    displaymatrix(sqrt1_Hii, "sqrt1_Hii");
  //    displaymatrix(wHalf, "wHalf");
}
