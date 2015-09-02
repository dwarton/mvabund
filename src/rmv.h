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
#ifndef _RMV_H
#define _RMV_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
/* GSL - GNU Scientific Library  */
/* #define GSL_CHECK_RANGE_OFF   */
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
/* ----------------------------- */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
/* ----------------------------- */
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
/* ------------------------------ */
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>

// MVN with positive-semi definite covariance matrix
int semirmvnorm(const gsl_rng *rnd, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result);

//multivariate normal distribution random number generator
int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result);
// multivariate normal density function
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);
// multivariate Student t distribution random number generator 
int rmvt(const gsl_rng *r, const int n, const gsl_vector *location, const gsl_matrix *scale, const int dof, gsl_vector *result);
// multivariate Student t density function
double dmvt(const int n, const gsl_vector *x, const gsl_vector *locationn, const gsl_matrix *scale, const int dof);
// Wishart distribution random number generator 
int rwishart(const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result);

#endif
