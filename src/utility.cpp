#include "utility.h"

int my_spmatrix_scale_row_ccs_double(gsl_spmatrix * m, gsl_vector * x) {
  if (m->size1 != x->size) {
      GSL_ERROR("x vector length does not match matrix", GSL_EBADLEN);
  } else {
    double* Ad = m->data;
	  size_t* Ai = m->i;
    size_t i;
    for (i = 0; i < m->nz; ++i) {
       double y = gsl_vector_get(x, Ai[i]);
       Ad[i] *= y;
    }
  }
  return GSL_SUCCESS;
}


