// utility functions
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// // 16-June-2011

#include "resampTest.h"

int testStatCalc(mv_mat *H0, mv_mat *H1, mv_Method *mmRef, const unsigned int ifcalcH1det, double *stat, gsl_vector *statj)
{
	unsigned int j;
        unsigned int nVars=H0->SS->size1;
        unsigned int nRows=H0->mat->size1;
	unsigned int nP0=H0->X->size2;
	unsigned int nP1=H1->X->size2;
	double scale = (double)(nRows-nP1)/(double)(nP1-nP0);
	double ss0j, ss1j, sum;
        double logDetSS0, logDetSS1;
	gsl_matrix *SS0 = gsl_matrix_alloc(nVars, nVars);
	gsl_matrix *SS1 = gsl_matrix_alloc(nVars, nVars);
	gsl_matrix *Sd = gsl_matrix_alloc(nVars, nVars);

        // univariate test
	for ( j=0; j<nVars; j++ ) {
            ss0j = gsl_matrix_get(H0->SS, j, j);
            ss1j = gsl_matrix_get(H1->SS, j, j);
	    if ( mmRef->test == LOGWILK ) 
	        gsl_vector_set(statj, j, nRows*(log(ss0j)-log(ss1j)));
            else // HOTELLING 
                gsl_vector_set(statj, j, (double)scale*(ss0j-ss1j)/ss1j);
	} 

	// multivariate test
        if ( mmRef->corr == IDENTITY ){ 
            sum = 0.0;
	    for (j=0; j<nVars;j++) sum += gsl_vector_get(statj, j);
        }
        else {
	    if (mmRef->corr == SHRINK) {
               gsl_matrix_memcpy(SS0, H0->R);
	       gsl_matrix_set_zero (Sd);
	       gsl_blas_dger (1.0, H0->sd, H0->sd, Sd);
	       gsl_matrix_mul_elements(SS0, Sd);
	       gsl_matrix_memcpy(SS1, H1->R);
	       gsl_matrix_set_zero (Sd);
	       gsl_blas_dger (1.0, H1->sd, H1->sd, Sd);
	       gsl_matrix_mul_elements(SS1, Sd);
            }
	    else {
	       gsl_matrix_memcpy(SS0, H0->SS);
	       gsl_matrix_memcpy(SS1, H1->SS);
	    }

            if (mmRef->test == LOGWILK) {
	      	logDetSS0 = logDet(SS0) ;		 
		 // used in summary to speed up 
	         if (ifcalcH1det == TRUE){
	            logDetSS1 = logDet(SS1);
		    H1->teststat = logDetSS1;
		 }
		 else logDetSS1 = H1->teststat;
                 sum = nRows * (logDetSS0 - logDetSS1);
            }
	    else {
	       // (SS0-SS1)SS1^-1
	       gsl_matrix_sub(SS0, SS1);
	       gsl_blas_dtrsm (CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, SS1, SS0);
	       sum = 0.0;
	       for (j=0; j<nVars;j++) sum += gsl_matrix_get(SS0, j, j); 
	       sum = sum*scale;
           }
	}
        // copy results
	*stat = sum;

	gsl_matrix_free(Sd);
	gsl_matrix_free(SS0);
	gsl_matrix_free(SS1);

	return 0;
}

int calcSS(gsl_matrix *Y, mv_mat *Hat, mv_Method *mmRef)
{
    gsl_set_error_handler_off();

    double dd;	
    unsigned int j;
    int status; 
    unsigned int nP=Hat->X->size2;
    unsigned int nRows=Hat->mat->size1;
    unsigned int nVars=Hat->SS->size1;
    gsl_matrix *tXX = gsl_matrix_alloc(nP, nP);
    gsl_matrix *tXXtX = gsl_matrix_alloc(nP, nRows);
    gsl_matrix *Sd = gsl_matrix_alloc(nVars, nVars);
     // It is possible later to feed data with more varialbes (columns) 
     // than observations (rows). So better use SVD than QR
     // As in manyglm use (X'X)^-1 X' and cholesky on (X'X) instead 
    
     gsl_matrix_set_identity(tXX);
     gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, Hat->X, 0.0, tXX);

     status=gsl_linalg_cholesky_decomp(tXX);
     if (status==GSL_EDOM) {
        gsl_matrix_set_identity(tXX);
        // add an eps*I to the diagonal matrix
        gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, Hat->X, TOL, tXX);
        status=gsl_linalg_cholesky_decomp(tXX);
     }
     gsl_linalg_cholesky_invert (tXX);
//     gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, Hat->X, Hat->X, 0.0, tXX); 
//     for (i=0; i<nRows; i++) {
//         gsl_vector_view hi=gsl_matrix_column(tXXtX, i);
//         gsl_vector_view xi=gsl_matrix_row(Hat->X, i);
//         GetCholstat(tXX, &xi.vector, &hi.vector);    
//     }

     gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, tXX, Hat->X, 0.0, tXXtX); 
     gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, tXXtX, Y, 0.0, Hat->Coef); 
     // Hat = X (XX)^-1 t(X) 
     gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Hat->X, tXXtX, 0.0, Hat->mat); 
     // Res = Y - Hat*Y 
     gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -1.0, Hat->mat, Y, 0.0, Hat->Res);  
     gsl_matrix_add (Hat->Res, Y);
     // subtractMean(Hat->Res);
     gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Hat->Res,Hat->Res,0.0,Hat->SS);  
     for (j=0; j<nVars; j++){
         dd = gsl_matrix_get(Hat->SS, j, j);
	 if (dd<TOL) dd = 0.5*M_1_PI;
	 gsl_matrix_set(Hat->SS, j, j, dd);
	 gsl_vector_set(Hat->sd, j, sqrt(dd));
     }
     if (mmRef->corr==IDENTITY) gsl_matrix_set_identity(Hat->R);
     else {
         gsl_matrix_memcpy(Hat->R, Hat->SS);
	 gsl_matrix_set_zero(Sd);
	 gsl_blas_dger (1.0, Hat->sd, Hat->sd, Sd);
	 gsl_matrix_div_elements(Hat->R, Sd);
	 if (mmRef->corr==SHRINK) {
	    // lambda * R + (1-lambda)*I
	    gsl_matrix_scale(Hat->R, (double)mmRef->shrink_param);
	    gsl_vector_view dj=gsl_matrix_diagonal(Hat->R);
	    gsl_vector_add_constant(&dj.vector, (double)(1.0-mmRef->shrink_param));
	  }
     }
 
     gsl_matrix_free(Sd);
     gsl_matrix_free(tXX);
     gsl_matrix_free(tXXtX);

     return 0;
}

double logDet(gsl_matrix *SS)
{
     // SS is a nVars x nVars real symmetric matric
     unsigned int nVars = SS->size1;

     // det(A) = prod(eig(A))
     gsl_eigen_symm_workspace *ws = gsl_eigen_symm_alloc(nVars);
     gsl_vector *eval=gsl_vector_alloc(nVars);    
     unsigned int j;
     double logDetSS=0.0;
     gsl_eigen_symm(SS, eval, ws);
 //    displayvector(eval, "SS eigen values");
     for ( j=0; j<nVars; j++ )
         logDetSS += log(gsl_vector_get(eval, j));

     gsl_eigen_symm_free(ws);
     gsl_vector_free(eval);
     return logDetSS;
  
/*
   // OR: det(A) = diag(LU)
     // fill SS
     gsl_matrix *LU=gsl_matrix_alloc(nVars, nVars);
     gsl_matrix_memcpy(LU, SS); 
     unsigned int i, j;
     int s;
//     for (i=0;i<nVars; i++)
//     for (j=i+1; j<nVars; j++)
//         gsl_matrix_set(LU, i, j, gsl_matrix_get(LU, j, i));

     gsl_permutation *p = gsl_permutation_alloc(nVars);
     gsl_linalg_LU_decomp (LU, p, &s);
     result = gsl_linalg_LU_det (LU, s);

     gsl_permutation_free(p);
     gsl_matrix_free(LU);

     return result;
*/
}

int is_sym_matrix(const gsl_matrix *mat)
{
     unsigned int i, j;
     if ( mat->size1 == mat->size2 ) {
        // test the upper triangle
        for ( i=0; i<mat->size1; i++)
            for ( j=i+1; j<mat->size2; j++ )
                if ( gsl_matrix_get(mat, i, j) != 0 )
	           return 0;
        return 1;
     }
     else
        return 0;
}

int subX(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    unsigned int j, k=0;
    unsigned int nP=ref->size;

    for (j=0; j<nP; j++) {
        if ( (unsigned int)gsl_vector_get(ref, j) > 0 ){
           gsl_vector_view col = gsl_matrix_column (X, j);
           gsl_matrix_set_col(Xi, k, &col.vector);
            k++;
	}    
     }
     return 0;
}

int subX1(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    unsigned int j, k=0;
    unsigned int nP=ref->size;

    for (j=0; j<nP; j++) {
        if ( (unsigned int)gsl_vector_get(ref, j) == 0 ){
            gsl_vector_view col = gsl_matrix_column (X, j);
            gsl_matrix_set_col(Xi, k, &col.vector);
            k++;
        }
    }

    return 0;
}


int subX2(gsl_matrix *X, unsigned int id, gsl_matrix *Xi)
{   
    unsigned int j, k=0;
    unsigned int nParam=X->size2;

    for (j=0; j<nParam; j++) {
        if ( j!=id ){
           gsl_vector_view col = gsl_matrix_column (X, j);
           gsl_matrix_set_col(Xi, k, &col.vector);
           k++;
	}   
    }
    return SUCCESS;
}

int subXrow(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    unsigned int j, k=0;
    unsigned int n=X->size1;

    for (j=0; j<n; j++) {       
        if ( (unsigned int)gsl_vector_get(ref, j)==0 ){
           gsl_vector_view row = gsl_matrix_row (X, j);
           gsl_matrix_set_row(Xi, k, &row.vector);
           k++;
        }
     }
     return 0;
}     

int subXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    unsigned int j, k=0;
    unsigned int n=X->size1;

    for (j=0; j<n; j++) {       
        if ( (unsigned int)gsl_vector_get(ref, j)>0 ){
           gsl_vector_view row = gsl_matrix_row (X, j);
           gsl_matrix_set_row(Xi, k, &row.vector);
           k++;
        }
     }
     return 0;
}     


int addXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    unsigned int j, k=0;
    unsigned int n=Xi->size1;

    gsl_matrix_set_zero(Xi);

    for (j=0; j<n; j++) {       
        if ( (unsigned int)gsl_vector_get(ref, j)>0 ){
           gsl_vector_view row = gsl_matrix_row (X, k);
           gsl_matrix_set_row(Xi, j, &row.vector);
           k++;
        }
     }
     return SUCCESS;
}     

int subXrow1(gsl_matrix *X, gsl_vector *ref0, gsl_vector *ref1, gsl_matrix *Xi)
{ 
    unsigned int j, k=0, id0, id1;
    unsigned int n=X->size1;

    for (j=0; j<n; j++) {       
        id0 = (unsigned int) gsl_vector_get(ref0, j);
        id1 = (unsigned int) gsl_vector_get(ref1, j);
        if ( (id0!=id1) & (id0==0) ){
           gsl_vector_view row = gsl_matrix_row (X, j);
           gsl_matrix_set_row(Xi, k, &row.vector);
           k++;
        }
     }
     return 0;
}     

int calcAdjustP(const unsigned int punit, const unsigned int nVars, double *bj, double *sj, double *pj, gsl_permutation *sortid)
{
    unsigned int k;
    
 //   printf("multiple test procedure [%d]\n", punit);
    if (punit == UNADJUST){
       for (k=0; k<nVars; k++)
       if (*(bj+k) >= *(sj+k))
          *(pj+k)=*(pj+k)+1;
     }
    else if (punit == SINGLESTEP) {
       gsl_vector_view bStatj=gsl_vector_view_array (bj, nVars);
       double maxstat=gsl_vector_max(&bStatj.vector);
       for (k=0; k<nVars; k++)
       if ( maxstat >= *(sj+k) )
           *(pj+k)=*(pj+k)+1;
     }
    else if (punit == FREESTEP){
       unsigned int sid, sid0=0;
       // successive maxima
       //gsl_permutation_fprintf(stdout, sortid, " %u");
 //      printf("\n");
       for (k=0; k<nVars; k++){
           sid = gsl_permutation_get(sortid, nVars-1-k);
//	   printf("%d ", (unsigned int)sid);
           if ( k>0 ) {
//	      printf("%d ", (unsigned int)sid0);
	      *(bj+sid)=MAX(*(bj+sid), *(bj+sid0));
           }   
//	   printf("%d ", (unsigned int)sid);
	   if ( *(bj+sid) >= *(sj+sid) ) 
	      *(pj+sid)=*(pj+sid)+1;
//	   if ( k> 0)    
//              printf("%.2f ", *(bj+sid0));   
//	   printf("%.2f ", *(bj+sid));   
//           printf("%.2f ", *(sj+sid));	  
//           printf("%d\n", (unsigned int)*(pj+sid));
	   sid0 = sid;   
    }  }
/*    else if (punit == STEPUP) {
       // successive minima
       for (k=0; k<nVars; k++){
           sid = gsl_permutation_get(sortid, k);
	   if ( k>0 ) {
	      *(bj+sid)=MIN(*(bj+sid), *(bj+sid0));
	   }
//           printf("%.3f ", (double)*(bj+sid));
	   //printf("%d ", (unsigned int)sid);
	   if (*(bj+sid) >= *(sj+sid))
	      *(pj+sid)=*(pj+sid)+1;
	   sid0 = sid;   
    }  }
    else if (punit == NOMONO) {
    }
//    printf("\n");
*/
    return 0;

}


int reinforceP(double *puj, unsigned int nVars, gsl_permutation *sortid)
{
// re-enforce monotonicity for univaraite tests
    unsigned int j, sid, sid0;
    for (j=1; j<nVars; j++) {
        sid=gsl_permutation_get(sortid, j);
	sid0=gsl_permutation_get(sortid, j-1);
	*(puj+sid) = MAX (*(puj+sid), *(puj+sid0));
    }	
    return 0;
}


int rcalc( gsl_matrix *Res, double shrink_param, unsigned int corr, gsl_matrix *SS)
{
    unsigned int j;
    unsigned int nRows = Res->size1;
    unsigned int nVars = Res->size2;
    
    gsl_vector *e = gsl_vector_alloc(nRows);
    gsl_vector_set_all (e, 1.0); 

    gsl_matrix_set_zero(SS);
    double mean, dd;	
    // subtract mean = Rs^T e(12) from res
    for ( j=0; j<nVars; j++ ) {
        gsl_vector_view rj=gsl_matrix_column(Res, j);
        gsl_blas_ddot (&rj.vector, e, &mean); 
        gsl_vector_add_constant (&rj.vector, -mean/nRows); 
        if ( corr == IDENTITY ) {
            gsl_blas_ddot (&rj.vector, &rj.vector, &dd); 
            if ( dd < 1e-10 ) dd = 0.5*M_1_PI;
            gsl_matrix_set(SS, j, j, dd);
        }
    } 
   // compute the covariance matrix SS= (res'*res)/Rows-1
    if ( corr != IDENTITY ){
        gsl_blas_dsyrk(CblasLower, CblasTrans,1.0, Res, 0.0, SS);
        gsl_matrix_scale ( SS, 1.0/(double)(nRows-1) );
        if ( corr == SHRINK ) {
           gsl_vector_view dj=gsl_matrix_diagonal (SS);
	   for (j=0; j<nVars; j++) {
	       dd = gsl_vector_get (&dj.vector, j);	       
	       if ( dd < 1e-10 ) // account for zero variances
	          gsl_vector_set (&dj.vector, j, 1.0/shrink_param); 
	       else  
	          gsl_vector_set (&dj.vector, j, dd/shrink_param); 
	   }    
        }
    } 
 
    gsl_vector_free(e);
    return 0;
}

int getHat(gsl_matrix *X, gsl_matrix *W, gsl_matrix *Hat) 
{
     double hii;
 
     unsigned int i, j;
     unsigned int nRows = X->size1;
     unsigned int nParam = X->size2;
     unsigned int nVars = Hat->size2;

     gsl_matrix *Xw = gsl_matrix_alloc(nRows, nParam);
     gsl_matrix *dW2 = gsl_matrix_alloc(nRows, nRows);
     gsl_matrix_set_zero(dW2);
     gsl_vector *t = gsl_vector_alloc(MIN(nRows, nParam));
     gsl_vector *twj = gsl_vector_alloc(nParam);
     gsl_vector *rj = gsl_vector_alloc(nRows);
     gsl_vector_view diagW, wj, qwj, xi;

     for ( j=0; j<nVars; j++) {
          gsl_matrix_memcpy (Xw, X);
          diagW = gsl_matrix_diagonal(dW2);
          wj = gsl_matrix_column(W, j);
          gsl_vector_memcpy(&diagW.vector, &wj.vector);
	  gsl_blas_dtrmm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, dW2, Xw); 
          gsl_linalg_QR_decomp (Xw, t); 
          // twj=(X'WX)^-1*XW = R\Q'wHalf where Xw = Q*R
          for (i=0; i<nRows; i++){
              qwj=gsl_matrix_column(dW2, i);
              gsl_linalg_QR_lssolve (Xw, t, &qwj.vector, twj, rj);

              xi=gsl_matrix_row(X, i);
	      gsl_blas_ddot(&xi.vector, twj, &hii);  //hii = xi*twj
              gsl_matrix_set(Hat, i, j, sqrt(1-hii));
          }         
     }
    
     gsl_matrix_free(Xw);
     gsl_matrix_free(dW2);
     gsl_vector_free(t);
     gsl_vector_free(twj);
     gsl_vector_free(rj);

     return 0;
}

int invLSQ(gsl_matrix *A, gsl_vector *b, gsl_vector *x)
{ 
   unsigned int dim1 = A->size1;
   unsigned int dim2 = A->size2;
   gsl_vector *t = gsl_vector_alloc(MIN(dim1, dim2));
   gsl_vector *r = gsl_vector_alloc(dim1);

   gsl_linalg_QR_decomp (A, t);
   gsl_linalg_QR_lssolve (A, t, b, x, r);

   gsl_vector_free(t);
   gsl_vector_free(r);
   
   return 0;
}


double GetSVDstat(gsl_matrix *A, gsl_vector *b, gsl_vector *work)
{ 
   unsigned int dim = A->size2;
   unsigned int k;
   double result;

   gsl_matrix *nullspace = gsl_matrix_alloc(dim, dim);
   gsl_vector *eigval = gsl_vector_alloc(dim);

   gsl_linalg_SV_decomp (A, nullspace, eigval, work);
   for (k=0; k<dim; k++) {
        if (gsl_vector_get(eigval, k) < TOL )
	    gsl_vector_set(eigval, k, 0);
   }	    
   gsl_linalg_SV_solve (A, nullspace, eigval, b, work);
   gsl_blas_ddot(b, work, &result);

   gsl_vector_free(eigval);
   gsl_matrix_free(nullspace);

   return result;
}


double GetCholstat(gsl_matrix *A, gsl_vector *b, gsl_vector *work)
{ 
   gsl_set_error_handler_off();

   unsigned int dim = A->size1;
   int status;
   double result;
   gsl_matrix *V = gsl_matrix_alloc(dim, dim);

   gsl_matrix_memcpy(V, A);
   status=gsl_linalg_cholesky_decomp(V);
   if (status==GSL_EDOM) {
      gsl_matrix_memcpy(V, A);
      gsl_vector_view diag_V = gsl_matrix_diagonal(V);
      gsl_vector_add_constant(&diag_V.vector, TOL);
      gsl_linalg_cholesky_decomp(V);
   }
   // result b^T V^-1 b
   gsl_linalg_cholesky_solve(V, b, work);
   gsl_blas_ddot(b, work, &result);
 
   gsl_matrix_free(V);

   return result;
}

// the manyglm version of rcalc
int GetR(gsl_matrix *Res, unsigned int corr, double lambda, gsl_matrix *R)
{
    if (corr == IDENTITY) {
       gsl_matrix_set_identity (R);
    }
    else {
       unsigned int j;
       unsigned int nRows = Res->size1;
       unsigned int nVars = Res->size2;
       double std;
       gsl_matrix *Sd =gsl_matrix_alloc(nVars, nVars);
       gsl_vector *s = gsl_vector_alloc(nVars);

       // get covariance matrix SS
       // note that residuals have already had means subtracted
       gsl_matrix_set_zero (R);
       gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0, Res, Res,0.0, R);
       gsl_matrix_scale (R, 1.0/(double)(nRows-1)); // SS

       gsl_matrix_set_all (Sd, 1.0);
       for ( j=0; j<nVars; j++ ) {               
           std = sqrt(gsl_matrix_get(R, j, j));          
           // account for zero variance
           gsl_vector_set(s, j, (std>TOL)?std:1/(2*M_1_PI)); 
       }
              
       gsl_blas_dsyr (CblasLower, 1.0, s, Sd);       
       gsl_matrix_set_zero(Sd);       
       gsl_blas_dger (1.0, s, s, Sd);
                     
       gsl_matrix_div_elements (R, Sd); // RR = SS./Sd              
       gsl_vector_view d = gsl_matrix_diagonal (R);       
       gsl_vector_set_all (&d.vector, 1.0); // ensure diag(R) = 1       
       //  displaymatrix(R, "corr-coefficents");
                                                        
       if ( corr == SHRINK ) { 
          gsl_matrix_scale (R, lambda);
          gsl_vector_add_constant(&d.vector,(double)(1.0-lambda));          
       }
                                                        
       gsl_matrix_free(Sd);       
       gsl_vector_free(s);       
       
    }
    
    return SUCCESS;
    
}

int subtractMean(gsl_matrix *dat)
{

    unsigned int nRows = dat->size1;
    unsigned int nVars = dat->size2;

//    gsl_matrix *t1, *Mean;
//    t1 = gsl_matrix_alloc(nRows, 1);
//    gsl_matrix_set_all (t1, 1.0);
//    Mean = gsl_matrix_alloc(nRows, nVars);
//    GetMean(t1, dat, Mean);
//    gsl_matrix_sub(dat, Mean);

    gsl_matrix *tmp = gsl_matrix_alloc(nRows, nVars);
    gsl_vector *t1 = gsl_vector_alloc(nRows);
    gsl_vector_set_all (t1, 1.0);
    
    gsl_vector *mean = gsl_vector_alloc(nVars);    
    gsl_blas_dgemv (CblasTrans, 1.0, dat, t1, 0.0, mean);
    gsl_vector_scale ( mean, (double)(1/nRows) );

    gsl_matrix_set_zero (tmp);
    gsl_blas_dger (1.0, t1, mean, tmp); 

    // subtract mean from dat
    gsl_matrix_sub(dat, tmp);

    gsl_vector_free(mean);
    gsl_vector_free(t1); 
    gsl_matrix_free(tmp);

    return SUCCESS;
}


