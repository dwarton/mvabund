#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "multi_thread.h"
#include <pthread.h>
#include "resampTest.h"
int resampAnovaCase(glm *model, gsl_matrix *bT, gsl_matrix *bX, gsl_matrix *bO, unsigned int i, mv_Method* tm, gsl_matrix* bootID, gsl_rng* rnd); 
int resampNonCase(glm *model, gsl_matrix *bT, unsigned int i, mv_Method* tm, gsl_matrix* bootID, gsl_rng* rnd, gsl_matrix* Xbeta, gsl_matrix* Sigma);
int GeeWald(glm *, gsl_matrix *, gsl_vector *, gsl_matrix *, mv_Method* tm);
int GeeScore(gsl_matrix *, glm *, gsl_vector *, gsl_matrix *, mv_Method* tm);
int GeeLR(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat, mv_Method* tm);

pthread_mutex_t anova_mutex = PTHREAD_MUTEX_INITIALIZER;
/*
int run_task(int total, int num_cores, void* task(void* data), void* data) {
   int ret;
   int i;
   pthread_t* threads = (pthread_t*) malloc(total* sizeof(pthread_t));
   pthread_attr_t attr;
   void *status;
   int task_num = 0;
   for( i = 0; i < total; i+= num_cores) {
      task_num = std::min(num_cores, total-num_cores);
      for (int j = 0; j < task_num; ++j) {
        thread_data* thd = new thread_data[1];
        thd->thread_id = (pthread_t) (i+j);
        thd->data  = data;
        ret = pthread_create(&threads[i+j], &attr, task, thd);
        if (ret) {
           printf("Error:unable to create thread %d err code %d\n", i+j, ret);
           return -1;
        }
      }
     for(int j = 0; j < task_num; j++) {
        ret = pthread_join(threads[i+j], &status);
        if (ret) {
           printf("Error:unable to join thread %d\n", i);
           return -1;
        }
     }
   }

  free(threads);
  return 0;
}

void *run_anova_mt(void *anova_pack) {
  anova_boot* data = (anova_boot*) anova_pack;
  GlmTest* glmtest =  data->gtest;
  glm* fit = data->fit;
  gsl_matrix* matrix = data->isXvarIn;
  glmtest->anova(fit, matrix);
  return NULL;
}
*/
void *anovaboot_mt(void *anova_pack) {
  anovaboot* th_data = (anovaboot*) anova_pack;
  glm* fit = th_data->fit;
  glm* PtrAlt = th_data->PtrAlt;
  glm* PtrNull = th_data->PtrNull;
  gsl_vector* ref0 = th_data->ref0;
  gsl_vector*  ref1 = th_data->ref1;
  gsl_matrix*  L1 = th_data->L1;
  int  ID0 = th_data->ID0;
  int ID1 = th_data->ID1;
  mv_Method* tm = th_data->tm;
  glm* bAlt = th_data->bAlt;
  glm* bNull = th_data->bNull;
  gsl_matrix* bY = th_data->bY;
  gsl_matrix* X1 = th_data->X1;
  gsl_matrix* X0 = th_data->X0;
  gsl_matrix*  bO = th_data->bO;
  gsl_matrix* BetaO = th_data->BetaO;
  gsl_matrix*  bootStore = th_data->bootStore;
	gsl_matrix* bootID = th_data->bootID;
	gsl_rng* rnd = th_data->rnd;
	gsl_matrix* XBeta = th_data->XBeta;
	gsl_matrix* Sigma = th_data->Sigma;

  int last = th_data->start_counter;
  int loop = th_data->loop_cnt;
    
  unsigned int nRows = tm->nRows, nVars = tm->nVars, nParam = tm->nParam;
  gsl_vector *bStat = gsl_vector_alloc(nVars + 1);
  gsl_matrix *Rlambda = gsl_matrix_alloc(tm->nVars, tm->nVars);
    
  for (int j = 0; j < loop; ++j) {
    gsl_vector_set_zero(bStat);
    if (tm->resamp == CASEBOOT) {
      resampAnovaCase(PtrAlt, bY, X1, bO, j, tm, bootID, rnd);
      subX(X1, ref0, X0);
    } else {
      resampNonCase(PtrNull, bY, j, tm, bootID, rnd, XBeta, Sigma);
      gsl_matrix_memcpy(bO, fit->Oref);
    }

    if (tm->test == WALD) {
      if (tm->resamp == CASEBOOT) {
        bAlt->regression(bY, X1, bO, NULL);
      } else {
        bAlt->regression(bY, NULL, bO, NULL);
      }
      double lambda = gsl_vector_get(tm->anova_lambda, ID1);
      GetR(bAlt->Res, tm->corr, lambda, Rlambda);
      GeeWald(bAlt, L1, bStat, Rlambda, tm);
    } else if (tm->test == SCORE) {
      if (tm->resamp == CASEBOOT) {
        bNull->regression(bY, X0, bO, NULL);
      } else {
        bNull->regression(bY, NULL, bO, NULL);
      }
      double lambda = gsl_vector_get(tm->anova_lambda, ID0);
      GetR(bNull->Res, tm->corr, lambda, Rlambda);
      GeeScore(X1, bNull, bStat, Rlambda, tm);
    } else {
      if (tm->resamp == CASEBOOT) {
        bNull->regression(bY, X0, bO, NULL);
      } else {
        bNull->regression(bY, NULL, bO, NULL);
      }
      addXrow2(bNull->Beta, ref1, BetaO);
      if (tm->resamp == CASEBOOT) {
        bAlt->regression(bY, X1, bO, BetaO);
      } else {
        bAlt->regression(bY, NULL, bO, BetaO);
      }
      GeeLR(bAlt, bNull, bStat, tm);
    }
		pthread_mutex_lock(&anova_mutex);
    gsl_matrix_set_row(bootStore, last + j, bStat);
		//std::cout << " boot last: "<< last  <<"  row num:" << last + j << std::endl;
 		pthread_mutex_unlock(&anova_mutex); 
  }
	gsl_matrix_free(Rlambda);
  gsl_vector_free(bStat);
  return NULL;
}


int resampAnovaCase(glm *model, gsl_matrix *bT, gsl_matrix *bX,
                             gsl_matrix *bO, unsigned int i, mv_Method* tm, gsl_matrix* bootID, gsl_rng* rnd) {
  gsl_set_error_handler_off();
  int status, isValid = TRUE;

  unsigned int j, id, nP;
  gsl_vector_view yj, xj, oj;
  nP = model->Xref->size2;
  gsl_matrix *tXX = gsl_matrix_alloc(nP, nP);
  unsigned int nRows = tm->nRows;

  while (isValid == TRUE) {
    for (j = 0; j < nRows; j++) {
      if (bootID != NULL)
        id = (unsigned int)gsl_matrix_get(bootID, i, j);
      else {
        if (tm->reprand == TRUE)
          id = (unsigned int)gsl_rng_uniform_int(rnd, nRows);
        else
          id = (unsigned int)nRows * Rf_runif(0, 1);
      }
      // resample Y and X and offset
      yj = gsl_matrix_row(model->Yref, id);
      xj = gsl_matrix_row(model->Xref, id);
      oj = gsl_matrix_row(model->Eta, id);
      // oj = gsl_matrix_row(model->Oref, id);
      gsl_matrix_set_row(bT, j, &yj.vector);
      gsl_matrix_set_row(bX, j, &xj.vector);
      gsl_matrix_set_row(bO, j, &oj.vector);
    }
    gsl_matrix_set_identity(tXX);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, bX, 0.0, tXX);
    status = gsl_linalg_cholesky_decomp(tXX);
    if (status != GSL_EDOM)
      break;
  }

  gsl_matrix_free(tXX);

  return SUCCESS;
}

int resampNonCase(glm *model, gsl_matrix *bT, unsigned int i, mv_Method* tm, gsl_matrix* bootID, gsl_rng* rnd, gsl_matrix* XBeta, gsl_matrix* Sigma) {
  unsigned int j, k, id;
  double bt, score, yij, mij;
  gsl_vector_view yj;
  unsigned int nRows = tm->nRows, nVars = tm->nVars;
  // to store Rf_unif
  //   gsl_vector *tmp = gsl_vector_alloc(nRows);
  //   gsl_permutation *vperm = gsl_permutation_alloc(nRows);
  double *tmp = (double *)malloc(nRows * sizeof(double));
  size_t* permid = NULL; 
	if (tm->resamp == PERMUTE) {
    permid = (size_t *)malloc(tm->nRows * sizeof(size_t));
    for (size_t i = 0; i < tm->nRows; i++)
      permid[i] = i;
  } else
    permid = NULL;

  // note that residuals have got means subtracted
  switch (tm->resamp) {
  case RESIBOOT:
    for (j = 0; j < tm->nRows; j++) {
      if (bootID != NULL)
        id = (unsigned int)gsl_matrix_get(bootID, i, j);
      else if (tm->reprand == TRUE)
        id = (unsigned int)gsl_rng_uniform_int(rnd, tm->nRows);
      else
        id = (unsigned int) tm->nRows * Rf_runif(0, 1);
      // bY = mu+(bootr*sqrt(variance))
      for (k = 0; k < nVars; k++) {
        bt = gsl_matrix_get(model->Mu, j, k) +
             sqrt(gsl_matrix_get(model->Var, j, k)) *
                 gsl_matrix_get(model->Res, id, k);
        bt = MAX(bt, 0.0);
        bt = MIN(bt, model->maxtol);
        gsl_matrix_set(bT, j, k, bt);
      }
    }
    break;
  case SCOREBOOT:
    for (j = 0; j < tm->nRows; j++) {
      if (bootID != NULL)
        score = (double)gsl_matrix_get(bootID, i, j);
      else if (tm->reprand == TRUE)
        score = gsl_ran_ugaussian(rnd);
      else
        score = Rf_rnorm(0.0, 1.0);
      // bY = mu + score*sqrt(variance)
      for (k = 0; k < nVars; k++) {
        bt = gsl_matrix_get(model->Mu, j, k) +
             sqrt(gsl_matrix_get(model->Var, j, k)) *
                 gsl_matrix_get(model->Res, j, k) * score;
        bt = MAX(bt, 0.0);
        bt = MIN(bt, model->maxtol);
        gsl_matrix_set(bT, j, k, bt);
      }
    }
    break;
  case PERMUTE:
    if (bootID == NULL) {
      if (tm->reprand == TRUE)
        gsl_ran_shuffle(rnd, permid, tm->nRows, sizeof(size_t));
      else { // Permutation with the randomness set in R
        for (j = 0; j < tm->nRows; j++)
          tmp[j] = Rf_runif(0, 1);
        gsl_sort_index(permid, tmp, 1, tm->nRows);
      }
    }
    for (j = 0; j < tm->nRows; j++) {
      if (bootID == NULL)
        id = permid[j];
      else
        id = (unsigned int)gsl_matrix_get(bootID, i, j);

      // bY = mu + bootr * sqrt(var)
      for (k = 0; k < nVars; k++) {
        bt = gsl_matrix_get(model->Mu, j, k) +
             sqrt(gsl_matrix_get(model->Var, j, k)) *
                 gsl_matrix_get(model->Res, id, k);
        bt = MAX(bt, 0.0);
        bt = MIN(bt, model->maxtol);
        gsl_matrix_set(bT, j, k, bt);
      }
    }
    break;
  case FREEPERM:
    if (bootID == NULL) {
      if (tm->reprand == TRUE)
        gsl_ran_shuffle(rnd, permid, nRows, sizeof(size_t));
      else { // Permutation with the randomness set in R
        for (j = 0; j < nRows; j++)
          tmp[j] = Rf_runif(0, 1);
        gsl_sort_index(permid, tmp, 1, nRows);
      }
    }
    for (j = 0; j < nRows; j++) {
      if (bootID == NULL)
        id = permid[j];
      else
        id = (unsigned int)gsl_matrix_get(bootID, i, j);

      yj = gsl_matrix_row(model->Yref, id);
      gsl_matrix_set_row(bT, j, &yj.vector);
    }
    break;
  case MONTECARLO:
    McSample(model, rnd, XBeta, Sigma, bT);
    break;
  case PITSBOOT:
    for (j = 0; j < nRows; j++) {
      if (bootID != NULL)
        id = (unsigned int)gsl_matrix_get(bootID, i, j);
      else if (tm->reprand == TRUE)
        id = (unsigned int)gsl_rng_uniform_int(rnd, nRows);
      else
        id = (unsigned int)Rf_runif(0, nRows);
      for (k = 0; k < nVars; k++) {
        bt = gsl_matrix_get(model->PitRes, id, k);
        mij = gsl_matrix_get(model->Mu, j, k);
        yij = model->cdfinv(bt, mij, model->theta[k]);
        gsl_matrix_set(bT, j, k, yij);
      }
    }
    break;
  default:
    GSL_ERROR("The resampling method is not supported", GSL_ERANGE);
    break;
  }

  free(tmp);
  if (permid != NULL)
    free(permid);

  return SUCCESS;
}


int GeeLR(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat, mv_Method* tm) {
  unsigned int nVars = tm->nVars;
  double val, result = 0;
  for (unsigned int j = 0; j < nVars; j++) { // univariates
    val = PtrAlt->ll[j] - PtrNull->ll[j];
    if (val < -0.1) {
      val = 0;
      if (tm->warning == TRUE)
        printf("Warning: Alt ll=%.4f < Null ll=%.4f\n", PtrAlt->ll[j],
               PtrNull->ll[j]);
    }
    gsl_vector_set(teststat, j + 1, val);
    result = result + val;
  }
  gsl_vector_set(teststat, 0, result); // multivariate
  return SUCCESS;
}

int GeeScore(gsl_matrix *X1, glm *PtrNull, gsl_vector *teststat, gsl_matrix* rlambda, mv_Method* tm) {
  gsl_set_error_handler_off();
	
	double eps = tm->tol;

  double result, alpha, sum = 0;
  unsigned int i, j, l, nP = X1->size2;
  unsigned int nVars = tm->nVars, nRows = tm->nRows;
  int status;

  gsl_vector *U = gsl_vector_alloc(nVars * nP);
  gsl_matrix *kRlNull = gsl_matrix_alloc(nVars * nP, nVars * nP);
  gsl_matrix_set_zero(kRlNull);
  gsl_matrix *XwX = gsl_matrix_alloc(nP, nP);
  gsl_vector *tmp = gsl_vector_alloc(nVars * nP);
  gsl_vector_view wj, uj, rj, tmp2; //, dj;
  gsl_matrix_view Rl;

  GrpMat *Z = (GrpMat *)malloc(nVars * sizeof(GrpMat));
  for (j = 0; j < nVars; j++) {
    Z[j].matrix = gsl_matrix_alloc(nRows, nP);
    // get W^1/2 * X
    wj = gsl_matrix_column(PtrNull->wHalf, j);
    for (i = 0; i < nP; i++)
      gsl_matrix_set_col(Z[j].matrix, i, &wj.vector);
    gsl_matrix_mul_elements(Z[j].matrix, X1);

    uj = gsl_vector_subvector(U, j * nP, nP);
    rj = gsl_matrix_column(PtrNull->Res, j);
    gsl_blas_dgemv(CblasTrans, 1, Z[j].matrix, &rj.vector, 0, &uj.vector);

    if ((tm->punit != NONE) || (tm->corr == IDENTITY)) {
      gsl_matrix_set_identity(XwX);
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1, Z[j].matrix, 0, XwX);
      status = gsl_linalg_cholesky_decomp(XwX);
      if (status == GSL_EDOM) {
        if (tm->warning == TRUE)
          printf("Warning: singular matrix in score test. An eps*I is added to "
                 "the singular matrix.\n");
        gsl_matrix_set_identity(XwX);
        gsl_blas_dsyrk(CblasLower, CblasTrans, 1, Z[j].matrix, eps, XwX);
        gsl_linalg_cholesky_decomp(XwX);
      }
      tmp2 = gsl_vector_subvector(tmp, 0, nP);
      gsl_linalg_cholesky_solve(XwX, &uj.vector, &tmp2.vector);
      gsl_blas_ddot(&uj.vector, &tmp2.vector, &result);
      gsl_vector_set(teststat, j + 1, result);
      sum = sum + result;
    }

    if (tm->corr != IDENTITY) {
      for (l = 0; l <= j; l++) { // lower half
        alpha = gsl_matrix_get(rlambda, j, l);
        Rl = gsl_matrix_submatrix(kRlNull, j * nP, l * nP, nP, nP);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, alpha, Z[j].matrix,
                       Z[l].matrix, 0, &Rl.matrix);
      }
    }
  } // end for j=1:nVars

  // multivariate test stat
  if (tm->corr == IDENTITY)
    gsl_vector_set(teststat, 0, sum);
  else {
    status = gsl_linalg_cholesky_decomp(kRlNull);
    if (status == GSL_EDOM) {
      if (tm->warning == TRUE)
        printf("Warning:singular kRlNull in multivariate score test.\n");
    }
    gsl_linalg_cholesky_solve(kRlNull, U, tmp);
    gsl_blas_ddot(U, tmp, &result);
    gsl_vector_set(teststat, 0, result);
  }

  // clear memory
  gsl_vector_free(U);
  gsl_vector_free(tmp);
  gsl_matrix_free(XwX);
  gsl_matrix_free(kRlNull);
  for (j = 0; j < nVars; j++)
    gsl_matrix_free(Z[j].matrix);
  free(Z);

  return SUCCESS;
}

// Wald Test used in both summary and anova (polymophism)
int GeeWald(glm *Alt, gsl_matrix *LL, gsl_vector *teststat, gsl_matrix* rlambda, mv_Method* tm) {
  gsl_set_error_handler_off();
	double eps = tm->tol;

  unsigned int i, j, l;
  double alpha, result, sum = 0;
  unsigned int nP = Alt->nParams;
  unsigned int nDF = LL->size1;
  unsigned int nVars = tm->nVars, nRows = tm->nRows;
  int status;
  gsl_vector *LBeta = gsl_vector_alloc(nVars * nDF);
  gsl_vector_set_zero(LBeta);
  gsl_matrix *w1jX1 = gsl_matrix_alloc(nRows, nP);
  gsl_matrix *XwX = gsl_matrix_alloc(nP, nP);
  gsl_matrix *Rl2 = gsl_matrix_alloc(nDF, nP);
  gsl_matrix *IinvN = gsl_matrix_alloc(nDF, nDF);
  gsl_matrix *IinvRl = gsl_matrix_alloc(nVars * nDF, nVars * nDF);
  gsl_vector *tmp = gsl_vector_alloc(nVars * nDF);
  gsl_vector_view tmp2, wj, LBj, bj; //, dj;
  gsl_matrix_view Rl;
  gsl_matrix_set_zero(IinvRl);
  GrpMat *Z = (GrpMat *)malloc(nVars * sizeof(GrpMat));
  for (j = 0; j < nVars; j++) {
    Z[j].matrix = gsl_matrix_alloc(nP, nRows);
    // w1jX1 = W^1/2 * X
    wj = gsl_matrix_column(Alt->wHalf, j);
    for (i = 0; i < nP; i++)
      gsl_matrix_set_col(w1jX1, i, &wj.vector);
    gsl_matrix_mul_elements(w1jX1, Alt->Xref);

    // LBeta = L*Beta
    LBj = gsl_vector_subvector(LBeta, j * nDF, nDF);
    bj = gsl_matrix_column(Alt->Beta, j);
    gsl_blas_dgemv(CblasNoTrans, 1, LL, &bj.vector, 0, &LBj.vector);

    // Z = (X^T W X)^-1 * X^T W^1/2.
    gsl_matrix_set_identity(XwX);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, w1jX1, 0.0, XwX);
    status = gsl_linalg_cholesky_decomp(XwX);
    if (status == GSL_EDOM) {
      if (tm->warning == TRUE)
        printf("Warning:singular matrix in wald test. An eps*I is added to the "
               "singular matrix.\n");
      gsl_matrix_set_identity(XwX);
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, w1jX1, eps, XwX);
      gsl_linalg_cholesky_decomp(XwX);
    }
    gsl_linalg_cholesky_invert(XwX);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, XwX, w1jX1, 0.0, Z[j].matrix);

    gsl_matrix_memcpy(Rl2, LL);
    gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, XwX,
                   Rl2); // L*(X'WX)^-1
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Rl2, LL, 0.0,
                   IinvN); // L*(X^T*W*X)^-1*L^T

    if ((tm->punit != NONE) || (tm->corr == IDENTITY)) {
      status = gsl_linalg_cholesky_decomp(IinvN);
      if (status == GSL_EDOM) {
        if (tm->warning == TRUE)
          printf("Warning:singular IinvN in wald test.\n");
      }
      tmp2 = gsl_vector_subvector(tmp, 0, nDF);
      gsl_linalg_cholesky_solve(IinvN, &LBj.vector, &tmp2.vector);
      gsl_blas_ddot(&LBj.vector, &tmp2.vector, &result);
      gsl_vector_set(teststat, j + 1, sqrt(result));
      sum = sum + result;
    }

    if (tm->corr != IDENTITY) {
      // IinvRl=L*vSandRl*L^T
      for (l = 0; l <= j; l++) {
        Rl = gsl_matrix_submatrix(IinvRl, j * nDF, l * nDF, nDF, nDF);
        alpha = gsl_matrix_get(rlambda, j, l);
        // borrow XwX space to store vSandRl
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, alpha, Z[j].matrix,
                       Z[l].matrix, 0.0, XwX);
        // Rl2 = L*vSandRl*L^T
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, LL, XwX, 0.0, Rl2);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Rl2, LL, 0.0, &Rl.matrix);
      } // end l
    }   // end if (tm->corr)

  } // end for j=1:nVars

  if (tm->corr == IDENTITY) {
    gsl_vector_set(teststat, 0, sqrt(sum));
	} else {
    status = gsl_linalg_cholesky_decomp(IinvRl);
    if (status == GSL_EDOM) {
      if (tm->warning == TRUE)
        printf("Warning:singular matrix in multivariate wald test.\n");
    }
    gsl_linalg_cholesky_solve(IinvRl, LBeta, tmp);
    gsl_blas_ddot(LBeta, tmp, &result);
    gsl_vector_set(teststat, 0, sqrt(result));
  }

  // free memory
  for (j = 0; j < nVars; j++)
    gsl_matrix_free(Z[j].matrix);
  free(Z);
  gsl_vector_free(LBeta);
  gsl_matrix_free(w1jX1);
  gsl_matrix_free(XwX);
  gsl_matrix_free(Rl2);
  gsl_matrix_free(IinvN);
  gsl_matrix_free(IinvRl);
  gsl_vector_free(tmp);

  return SUCCESS;
}


