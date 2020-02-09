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
#include "resampTest.h"
/*
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

struct anova_boot {
  GlmTest* gtest;
	glm *fit;
	gsl_matrix *isXvarIn;
};
void *run_anova_mt(void *anova_pack);

int run_task(int total, int num_cores, void* task(void* data), void* data);
*/
struct anovaboot {

//read only vars
glm* fit;
glm* PtrAlt;
glm* PtrNull;
gsl_vector* ref0;
gsl_vector* ref1;
gsl_matrix* L1;
int ID0;
int ID1;
mv_Method* tm;
//glmtest vars

gsl_rng *rnd;
gsl_matrix *bootID;

gsl_matrix *XBeta;
gsl_matrix *Sigma; // used in monte carlo simulation

//R/W vars
glm* bAlt;
glm* bNull;
gsl_matrix* bY;
gsl_matrix* X1;
gsl_matrix* X0;
gsl_matrix* bO;
gsl_matrix* BetaO;


//mux loc vars
gsl_matrix* bootStore;

int start_counter;
int loop_cnt;
};

void *anovaboot_mt(void *anova_pack);

#endif

