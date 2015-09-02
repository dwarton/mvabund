// more utility functions, mostly used in glmtest or simulations
// Author: Yi Wang (yi dot want at unsw dot edu dot au)
// 20-Apr-2011

#include "resampTest.h"

int GetMean(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *Mu)
{
    unsigned int nRows = X->size1;
    unsigned int nParam = X->size2;
    unsigned int nVars = Y->size2;
    unsigned int j;

    gsl_matrix *Coef=gsl_matrix_alloc(nParam, nVars);
    gsl_matrix *U=gsl_matrix_alloc(nRows, nParam);
    gsl_vector *t=gsl_vector_alloc(MIN(nRows, nParam));
    
    gsl_matrix_memcpy(U, X);
    gsl_linalg_QR_decomp (U, t);
    gsl_vector_view yj, cj, rj;
    // Compute coef and residual from Y and QR(X)
    for (j=0; j<nVars; j++){
        yj = gsl_matrix_column(Y, j);
        cj = gsl_matrix_column(Coef, j);
        rj = gsl_matrix_column(Mu, j);
	gsl_linalg_QR_lssolve (U, t, &yj.vector, &cj.vector, &rj.vector);
    }
    // Mu=Y-residual
    gsl_matrix_sub (Mu, Y);
    gsl_matrix_scale (Mu, -1.0);

    gsl_matrix_free(Coef);
    gsl_matrix_free(U);
    gsl_vector_free(t);

    return 0;
}

int GetCov (gsl_matrix *Mu, gsl_matrix *Y, unsigned int AR1MAT, gsl_matrix *Sigma)
{
    unsigned int nRows=Y->size1;
    unsigned int nVars=Y->size2;
    unsigned int i, j;
    double tmp, d1, d2;

    /* ---------------------------- */
    /*   Sample covariance matrix   */
    /* ---------------------------- */
    gsl_matrix *Res=gsl_matrix_alloc (nRows, nVars);
    gsl_matrix *SS=gsl_matrix_alloc (nVars, nVars);
    gsl_matrix_memcpy (Res, Y);
    gsl_matrix_sub (Res, Mu);
    gsl_matrix_set_zero (SS);
    // SS = RES^T RES / (nRows-1)
    gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, Res, 0.0, SS);
/*    for (i=0; i<nVars; i++) 
    for (j=0; j<nVars; j++) {
        // copy lower triangle to upper triangle
        gsl_matrix_set(SS, i, j, gsl_matrix_get(SS, j, i));
    }    
    // apply shrinkage here ?
*/
    /* ---------------------------- */
    /*      Covariance matrix       */
    /* ---------------------------- */
    gsl_matrix_set_identity (Sigma);  
    gsl_vector_view dS = gsl_matrix_diagonal (Sigma);
    gsl_vector_view sig = gsl_matrix_diagonal (SS);
    if (AR1MAT == 0) {  
       // Sigma=sigma*I       
       gsl_vector_memcpy (&dS.vector, &sig.vector);
       gsl_vector_scale (&dS.vector, (double)1.0/(nRows-1));
    }
    else if (AR1MAT == 1) { 
       // Sigma = dd^1/2*AR1*dd^1/2
       for (i=0; i<nVars; i++) {
	   // AR1 structure
           for (j=i+1; j<nVars; j++){
	       tmp = gsl_matrix_get(Sigma, i, j-1);
               gsl_matrix_set(Sigma, i, j, tmp*RHO);
	       gsl_matrix_set(Sigma, j, i, tmp*RHO);
           } 
       }
       // dd 
       // copy variances from SS
       gsl_vector_memcpy (&dS.vector, &sig.vector);
       for (i=0; i<nVars; i++) {
          d1=sqrt(gsl_matrix_get(SS, i, i)); 
	  for (j=i+1; j<nVars; j++){
	      d2=sqrt(gsl_matrix_get(SS, j, j));
	      tmp = d1*d2/(nRows-1); // di*dj
	      gsl_matrix_set(Sigma, i, j, gsl_matrix_get(Sigma, i, j)*tmp);
              gsl_matrix_set(Sigma, j, i, gsl_matrix_get(Sigma, j, i)*tmp);
//	      printf("%.2f ", tmp);
          }
       }
    }
    else if (AR1MAT == 2) {
        // Sigma = (n-1)/n S | H0 )
        gsl_matrix_memcpy (Sigma, SS);
        gsl_matrix_scale (Sigma, (double)1.0/nRows);
    }

    gsl_matrix_free (Res);
    gsl_matrix_free (SS);

    return SUCCESS;
}

int GetMeanCov(gsl_matrix *X, gsl_matrix *Y, mv_Method *mm, unsigned int AR1MAT, gsl_matrix *Mu, gsl_matrix *Sigma)
{
    unsigned int i, j;
    double tmp, d1, d2;
    unsigned int nRows = X->size1; 
    unsigned int nParam = X->size2;
    unsigned int nVars = Y->size2;
    unsigned int INCOR = mm->corr;
    mv_mat H;
    H.SS = gsl_matrix_alloc (nVars, nVars);
    H.mat = gsl_matrix_alloc(nRows, nRows);
//    H.Coef = gsl_matrix_alloc(nParam, nVars);
    H.Res = gsl_matrix_alloc(nRows, nVars);
    H.X = gsl_matrix_alloc(nRows, nParam);
    gsl_matrix_memcpy(H.X, X);
    // fit the model to Y
    mm->corr = NOSHRINK;  
                      // hat, coef, SS
    calcSS (Y, &H, mm); 
    mm->corr = INCOR;

    /* ---------------------- */
    /*   Mean matrix under H  */
    /* ---------------------- */
    gsl_matrix_memcpy (Mu, Y);
    // Mu = H.Y = Y - H.Res
    gsl_matrix_sub (Mu, H.Res);
//    displaymatrix (Mu, "Mu");

    /* ---------------------------- */
    /*   Covariance matrix under H  */
    /* ---------------------------- */
    gsl_vector_view sig = gsl_matrix_diagonal (H.SS);

    gsl_matrix_set_identity (Sigma);  
    gsl_vector_view dS = gsl_matrix_diagonal (Sigma);
    if (AR1MAT == 0) {  
       // Sigma=sigma*I       
       gsl_vector_memcpy (&dS.vector, &sig.vector);
    }
    else if (AR1MAT == 1) { 
       // Sigma = dd^1/2*AR1*dd^1/2
       for (i=0; i<nVars; i++) {
	   // AR1 structure
           for (j=i+1; j<nVars; j++){
	       tmp = gsl_matrix_get(Sigma, i, j-1);
               gsl_matrix_set(Sigma, i, j, tmp*RHO);
	       gsl_matrix_set(Sigma, j, i, tmp*RHO);
           } 
       }
       // dd 
       // copy variances from H0.SS
       gsl_vector_memcpy (&dS.vector, &sig.vector);
       for (i=0; i<nVars; i++) {
          d1=sqrt(gsl_matrix_get(H.SS, i, i)); 
	  for (j=i+1; j<nVars; j++){
	      d2=sqrt(gsl_matrix_get(H.SS, j, j));
	      tmp = d1*d2; // di*dj
	      gsl_matrix_set(Sigma, i, j, gsl_matrix_get(Sigma, i, j)*tmp);
              gsl_matrix_set(Sigma, j, i, gsl_matrix_get(Sigma, j, i)*tmp);
//	      printf("%.2f ", tmp);
          }
//          printf("\n");	   
       }
    }
    else if (AR1MAT == 2) {
        // Sigma = (n-1)/n S | H0 )
        gsl_matrix_memcpy (Sigma, H.SS);
        gsl_matrix_scale (Sigma, (double)(nRows-1)/nRows);
	// copy lower triangle to upper triangle
        for (i=0;i<nVars; i++)
        for (j=i+1; j<nVars; j++)
            gsl_matrix_set(Sigma, i, j, gsl_matrix_get(Sigma, j, i));
    }
    else if (AR1MAT == 3) { 
        // sigma*I
    }
    else if (AR1MAT == 4) { // Sigma = AR1
       for (i=0; i<nVars; i++)  
       for (j=i+1; j<nVars; j++){
           tmp = gsl_matrix_get(Sigma, i, j-1);
           gsl_matrix_set(Sigma, i, j, tmp*RHO);
           gsl_matrix_set(Sigma, j, i, tmp*RHO);
       }
    }
    else if (AR1MAT == 5) { // Sigma = R 
        // dd
        gsl_vector_memcpy (&dS.vector, &sig.vector);
	for (i=0; i<nVars; i++) {
           d1=sqrt(gsl_matrix_get(H.SS, i, i)); 
           for (j=i+1; j<nVars; j++){
	       d2=sqrt(gsl_matrix_get(H.SS, j, j));  // consider chol
               tmp = d1*d2; // di*dj
               gsl_matrix_set(Sigma, i, j, tmp);
               gsl_matrix_set(Sigma, j, i, tmp);
	       // copy lower triangle to upper triangle
               gsl_matrix_set(H.SS, i, j, gsl_matrix_get(H.SS, j, i));
	   }
        }
	// R = SS./dd
	gsl_matrix_div_elements (H.SS, Sigma);
	gsl_matrix_memcpy(Sigma, H.SS);
    }
//    displaymatrix(Sigma, "Covariance Matrix");

    gsl_matrix_free(H.mat);
//    gsl_matrix_free(H.Coef);
    gsl_matrix_free(H.SS);
    gsl_matrix_free(H.X);
    gsl_matrix_free(H.Res);

    return 0;
}

int GetPdstbtion(double *p, unsigned int nVars, unsigned int *isH0var, unsigned int *cnt, unsigned int *cntfwe)
{
    // The multivariate test    
    if ( *p < ALFA+TOL) cnt[0]=cnt[0]+1;
//    printf("%.2f(%d) ", *p, cnt[0]);	
    // The univariate test
    double minP = 1.0;
    for ( unsigned int j=1; j<nVars+1; j++ ) {               
	if ( *(p+j) < ALFA+TOL ) cnt[j]=cnt[j]+1;	
//	printf("%.2f(%d) ", *(p+j), cnt[j]);	
        // strong FWE control over the null variables only	
	if ( (isH0var[j-1]==TRUE) & (*(p+j)<minP) ) minP=*(p+j);
    }    
    if ( minP < ALFA+TOL ) *cntfwe=*cntfwe+1;

//    printf("\n minP=%.2f, cntfwe=%d \n", minP, *cntfwe);
    return 0;
}   

// ----------------------------------------------------- //
//  Select test variables according to variances under H0 //
// ----------------------------------------------------- //
int GetH0var(gsl_matrix *Sigma, unsigned int *isH0var)
{
    unsigned int j, h0id, h1id; 
    unsigned int nVars = Sigma->size1;
    size_t *srtid = (size_t *)malloc(nVars*sizeof(size_t));

    // Sort variances in descending order
    gsl_vector_view sig=gsl_matrix_diagonal (Sigma);
    gsl_sort_vector_largest_index (srtid, nVars, &sig.vector);

    for (j=0; j<nVars; j+=2) {
        h0id = srtid[j];
        isH0var[h0id] = TRUE;
        h1id = srtid[j+1];
        isH0var[h1id] = FALSE;
    }      
/*   
    for (j=0; j<nVars/2; j++){
        h0id = srtid[j];
        isH0var[h0id] = TRUE;
    }     
    for (j=nVars/2; j<nVars; j++){
        h1id = srtid[j];
        isH0var[h1id] = FALSE;
    }     
*/
    free(srtid);
    return 0;
}

int setMonteCarlo(glm *model, gsl_matrix *XBeta, gsl_matrix *Sigma)
{
   unsigned int j;
   unsigned int mtype=model->mmRef->model;
   double k, vij, sd, scale;
 
   gsl_matrix_memcpy (XBeta, model->Eta);
   if (mtype == POISSON) {
       gsl_matrix_set_identity (Sigma);
       // Assuming no random effects, i.e. e*=0
   }
   else if (mtype == BIN) {
       if (model->n==1) {
       // Adjusting all betas and assuming var=1 
       // logit(M) = X * sqrt(1+0.346 var) beta = 1.1601 Eta 
       // See MATH5885 LDA lecture notes W9-11, Section 6.6.1
          sd = 1;
          k = 16*sqrt(double(3))/15/M_PI;
          scale = sqrt(1 + gsl_pow_2(k)*gsl_pow_2(sd));
          gsl_matrix_scale (XBeta, scale);
       }
   }
   else if (model->mmRef->model == NB) {
       // Adjusting the intercept to account for random effects
       //  i.e., M = X * beta - 0.5 * var 
       // var = log(1+phi)    
       gsl_matrix *Sd = gsl_matrix_alloc (model->nVars, model->nVars);
       gsl_vector  *s = gsl_vector_alloc (model->nVars);
       gsl_vector_view mj, d;

       for ( j=0; j<model->nVars; j++) {
           mj=gsl_matrix_column (XBeta, j);
           // adjust E(mj) = X*beta for the random effects
           vij = log(1+model->theta[j])-log(model->theta[j]);
           gsl_vector_add_constant(&mj.vector, -0.5*vij);
           gsl_vector_set(s, j, sqrt(vij));
       }
       gsl_matrix_set_zero (Sd);
       gsl_blas_dger (1.0, s, s, Sd);

       // if phi=0, then vij=0, i.e., no random effects (independence)
       // So it has zero impact / correlation on other variables 
       d = gsl_matrix_diagonal(Sd);
       for (j=0; j<model->nVars; j++) 
           if (model->theta[j]>100) gsl_vector_set(&d.vector, j, 1.0);

       // Sigma = diag(var)*R*diag(var)
       gsl_matrix_mul_elements(Sigma, Sd);
       // displaymatrix(Sigma, "log-normal Sigma");

       // free memory
       gsl_matrix_free(Sd);
       gsl_vector_free(s);
   }
   else 
       GSL_ERROR("The model type is not supported", GSL_ERANGE);

   return SUCCESS;
}

int McSample(glm *model, gsl_rng *rnd, gsl_matrix *XBeta, gsl_matrix *Sigma, gsl_matrix *bY) 
{
   unsigned int j, k;
   unsigned int nRows=XBeta->size1;
   unsigned int nVars=Sigma->size1;
   gsl_vector_view yj;
   double eij, mij, yij;

   if (model->mmRef->model == NB){  // Poisson log-normal
      for (j=0; j<nRows; j++) {
          yj = gsl_matrix_row(bY, j);
          semirmvnorm(rnd, nVars, Sigma, &yj.vector); // random effect
          for (k=0; k<nVars; k++) {
              eij=gsl_matrix_get (XBeta, j, k);
              // m_j = X_j * Beta_j + random_effect
              if ( model->theta[k]>100 ) { // add random effect
                 eij = eij + gsl_vector_get(&yj.vector, k);
              }
              // Sample from Poisson-Log-Normal dist
              yij = Rf_rpois(exp(eij));
              gsl_matrix_set(bY, j, k, yij);
          }
       }
   }
   else if (model->mmRef->model==BIN) {
       for (j=0; j<nRows; j++) {
            yj = gsl_matrix_row(bY, j);
            semirmvnorm(rnd, nVars, Sigma, &yj.vector); // random effect
            for (k=0; k<nVars; k++) {
                eij = gsl_matrix_get (XBeta, j, k); // logit(m)
                eij = eij + gsl_vector_get(&yj.vector, k);
                mij = model->invLink(eij);
                yij = model->genRandist(mij, model->theta[k]);
                gsl_matrix_set(bY, j, k, yij);
            }
       }
   }
   else {
       // Method 1 use R random gen func directly, no random effects
       for (j=0; j<nRows; j++)
       for (k=0; k<nVars; k++) {
            mij = gsl_matrix_get(model->Mu, j, k);
            yij = model->genRandist(mij, model->theta[k]);
            gsl_matrix_set(bY, j, k, yij);
       }
   }

   return SUCCESS;
}


