#ifndef RUNTHREAD
#define RUNTHREAD
#include <pthread.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>
#include <math.h>

struct thread_data {
   int  thread_id;
   void *data;
};
struct mt_data{
	  gsl_vector_view y;
		gsl_vector_view m;
		double* dl;
		double* ddl;
		double k;
};

int run_task(int total, int num_cores, void* task(void* data), void* data);


#endif

