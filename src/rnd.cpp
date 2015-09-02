/***************************************************************************************
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 *  Wishart random number generator
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  AUTHOR 
 *     Ralph dos Santos Silva,  [EMAIL PROTECTED]
 *     March, 2006       
***************************************************************************************/
//#include "rmv.h"

/***********************************************************************/
/* multivariate normal distribution random number generator 
   n	   dimension of the random vetor
   mean	   vector of means of size n
   var	   variance matrix of dimension n x n
   result  output variable with a sigle random vector normal distribution generation
************************************************************************/

#include "resampTest.h"

int rmvnorm(const gsl_rng *r, const unsigned int n, const gsl_matrix *Sigma, gsl_vector *randeffect)
{
    unsigned int k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);

    gsl_matrix_memcpy(work, Sigma);
    gsl_linalg_cholesky_decomp(work);

    for(k=0; k<n; k++)
	gsl_vector_set(randeffect, k, gsl_ran_ugaussian(r) );

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, randeffect);

    gsl_matrix_free(work);

    return 0;
}

/*********************************************************/
/* MVN with positive semi-definite covariance matrix     */
/* i.e. var is rank deficiant (has zero eigen values     */
/* use eigen decomposition instead of cholesky           */
/* see cholcov.m in matlab                               */
/*********************************************************/
int semirmvnorm(const gsl_rng *rnd, const unsigned int n, const gsl_matrix *Sigma, gsl_vector *randeffect)
{
    unsigned int k, r=0;
    double lambda;
    gsl_matrix *work = gsl_matrix_alloc(n,n);

    gsl_matrix_memcpy(work, Sigma);
//    replace cholesky with eigen decomposition
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
    gsl_vector *eval=gsl_vector_alloc (n);
    gsl_matrix *evec=gsl_matrix_alloc (n, n);
    // work = evec*diag(eval)*t(evec)
    gsl_eigen_symmv (work, eval, evec, w);
//    displayvector (eval, "eigen values of work");
//    displaymatrix (evec, "eigen vector of work");
    for (k=0; k<n; k++) {
        gsl_vector_view evec_i=gsl_matrix_column(evec, k);
        lambda=gsl_vector_get(eval, k);
        if (lambda>10e-10){ // non-zero variables
           // U = t(eval(r)*evec(:, r))
	   gsl_vector_scale (&evec_i.vector, sqrt(lambda));
	   // copy U to work 
	   gsl_matrix_set_col(work, r, &evec_i.vector);
	   r++;
	}   
    }
//    printf("r=%d.\n", r);
    gsl_matrix_view U=gsl_matrix_submatrix (work, 0, 0, n, r);
//    displaymatrix (&U.matrix, "partial eigen vectors");    

    // generate standard normal vector  
    gsl_vector *z=gsl_vector_alloc(r);
    for(k=0; k<r; k++)
	gsl_vector_set( z, k, gsl_ran_ugaussian(rnd) );
//    displayvector (z, "z"); 
    // X_i = mu_i + t(U)*z 
    gsl_blas_dgemv (CblasNoTrans, 1.0, &U.matrix, z, 0.0, randeffect);
//    displayvector (randeffect, "randeffect");

    gsl_matrix_free(work);
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_vector_free(z);

    return 0;
}


/*************************************************************************/
/* multivariate normal density function    */
/*
*	n	dimension of the random vetor
*	mean	vector of means of size n
*	var	variance matrix of dimension n x n
*/
/*************************************************************************/
double dmvnorm(const unsigned int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var)
{
   int s;
   double ax,ay;
   gsl_vector *ym, *xm;
   gsl_matrix *work = gsl_matrix_alloc(n,n), 
              *winv = gsl_matrix_alloc(n,n);
   gsl_permutation *p = gsl_permutation_alloc(n);

   gsl_matrix_memcpy( work, var );
   gsl_linalg_LU_decomp( work, p, &s );
   gsl_linalg_LU_invert( work, p, winv );
   ax = gsl_linalg_LU_det( work, s );
   gsl_matrix_free( work );
   gsl_permutation_free( p );

   xm = gsl_vector_alloc(n);
   gsl_vector_memcpy( xm, x);
   gsl_vector_sub( xm, mean );
   ym = gsl_vector_alloc(n);
   gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
   gsl_matrix_free( winv );
   gsl_blas_ddot( xm, ym, &ay);
   gsl_vector_free(xm);
   gsl_vector_free(ym);
   ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),double(n))*ax );

   return ay;
}

/*************************************************************************************/
/* multivariate Student t distribution random number generator */
/*
*	n	 dimension of the random vetor
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*	result	 output variable with a single random vector normal distribution generation
*/
/*************************************************************************************/
int rmvt(const gsl_rng *r, const unsigned int n, const gsl_vector *location, const gsl_matrix *scale, const unsigned int dof, gsl_vector *result)
{
    unsigned int k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);
    double ax = 0.5*dof; 

    ax = gsl_ran_gamma(r,ax,(1/ax));     /* gamma distribution */

    gsl_matrix_memcpy(work,scale);
    gsl_matrix_scale(work,(1/ax));       /* scaling the matrix */
    gsl_linalg_cholesky_decomp(work);

    for(k=0; k<n; k++)
	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
    gsl_vector_add(result, location);

    gsl_matrix_free(work);

    return 0;
}

/*****************************************************************************************************************/
/* multivariate Student t density function */
/*
*	n	 dimension of the random vetor
*	location vector of locations of size n
*	scale	 scale matrix of dimension n x n
*	dof	 degrees of freedom
*/
/*****************************************************************************************************************/
double dmvt(const unsigned int n, const gsl_vector *x, const gsl_vector *location, const gsl_matrix *scale, const unsigned int dof)
{
    int s;
    double ax,ay,az=0.5*(dof + n);
    gsl_vector *ym, *xm;
    gsl_matrix *work = gsl_matrix_alloc(n,n), 
               *winv = gsl_matrix_alloc(n,n);
    gsl_permutation *p = gsl_permutation_alloc(n);

    gsl_matrix_memcpy( work, scale );
    gsl_linalg_LU_decomp( work, p, &s );
    gsl_linalg_LU_invert( work, p, winv );
    ax = gsl_linalg_LU_det( work, s );
    gsl_matrix_free( work );
    gsl_permutation_free( p );

    xm = gsl_vector_alloc(n);
    gsl_vector_memcpy( xm, x);
    gsl_vector_sub( xm, location );
    ym = gsl_vector_alloc(n);
    gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
    gsl_matrix_free( winv );
    gsl_blas_ddot( xm, ym, &ay);
    gsl_vector_free(xm);
    gsl_vector_free(ym);

    ay = pow((1+ay/dof),-az)*gsl_sf_gamma(az)/(gsl_sf_gamma(0.5*dof)*sqrt( pow((dof*M_PI),double(n))*ax ));

    return ay;
}

/***************************************************************************************/
/* Wishart distribution random number generator */
/*
*	n	 gives the dimension of the random matrix
*	dof	 degrees of freedom
*	scale	 scale matrix of dimension n x n
*	result	 output variable with a single random matrix Wishart distribution generation
*/
/***************************************************************************************/
int rwishart(const gsl_rng *r, const unsigned int n, const unsigned int dof, const gsl_matrix *scale, gsl_matrix *result)
{
    unsigned int k,l;
    gsl_matrix *work = gsl_matrix_calloc(n,n);

    for(k=0; k<n; k++){
	gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
	for(l=0; l<k; l++)
	    gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
    }
    gsl_matrix_memcpy(result,scale);
    gsl_linalg_cholesky_decomp(result);
    gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,result,work);
    gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,result);

    return 0;
}
