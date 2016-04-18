// Resampling-based hypothesis test for comparing multivariate linear models
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 16-Nov-2009


#include "resampTest.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h> 
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort_double.h>

// Use gsl functions as much as poosible to increase speed and stabilitiy

AnovaTest::AnovaTest(mv_Method *mm, gsl_matrix *Y, gsl_matrix *X, gsl_matrix *isXvarIn):mmRef(mm), Yref(Y), Xref(X), inRef(isXvarIn)
{
    unsigned int hid, aid;
    unsigned int i, j, count;
    nModels=inRef->size1, nParam=Xref->size2;
    nRows=Yref->size1, nVars=Yref->size2; 

//  printf("initialize public variables: stats\n");
    multstat=(double *)malloc((nModels-1)*sizeof(double));
    Pmultstat = (double *)malloc((nModels-1)*sizeof(double));
    for (j=0; j<nModels-1; j++) *(Pmultstat+j)=0.0; 
    dfDiff = (unsigned int *)malloc((nModels-1)*sizeof(unsigned int));

    statj = gsl_matrix_alloc(nModels-1, nVars);
    Pstatj = gsl_matrix_alloc(nModels-1, nVars);
    gsl_matrix_set_zero(Pstatj);

    bStatj = gsl_vector_alloc(nVars);
    Hats = (mv_mat *)malloc(nModels*sizeof(mv_mat)); 
    sortid = (gsl_permutation **)malloc((nModels-1)*sizeof(gsl_permutation *));
    
    for (i=0; i<nModels; i++ ) {
        // Hats[i]
        Hats[i].mat=gsl_matrix_alloc(nRows, nRows);
        Hats[i].SS=gsl_matrix_alloc(nVars, nVars);
        Hats[i].R=gsl_matrix_alloc(nVars, nVars);
        Hats[i].Res=gsl_matrix_alloc(nRows, nVars);
        Hats[i].Y = gsl_matrix_alloc(nRows, nVars);
        Hats[i].sd = gsl_vector_alloc(nVars);
	count = 0;
	for (j=0; j<nParam; j++){
	    count+=(unsigned int)gsl_matrix_get(inRef, i, j);
	}
//	printf("count=%d \n", count);
	Hats[i].X = gsl_matrix_alloc(nRows, count);
	Hats[i].Coef=gsl_matrix_alloc(count, nVars);
        gsl_vector_view refi=gsl_matrix_row(inRef, i);
	subX(Xref, &refi.vector, Hats[i].X);
        calcSS(Yref, &(Hats[i]), mmRef);
//	displaymatrix(Hats[i].SS, "SS");
    }

    for (i=1; i<nModels; i++) {
        hid = i; aid = i-1;
        if ( mmRef->resamp != CASEBOOT ) {
            // fit = Y- resi 
            gsl_matrix_memcpy (Hats[i].Y, Yref);
            gsl_matrix_sub (Hats[i].Y, Hats[i].Res);
        } 
        gsl_vector_view statij = gsl_matrix_row(statj, aid);
        testStatCalc(&(Hats[hid]), &(Hats[aid]), mmRef, TRUE, (multstat+aid), &statij.vector); 
	dfDiff[aid] = Hats[aid].X->size2-Hats[hid].X->size2;
        // sortid
        sortid[aid] = gsl_permutation_alloc(nVars);
        gsl_sort_vector_index (sortid[aid], &statij.vector); 
        // rearrange sortid in descending order
        gsl_permutation_reverse (sortid[aid]);
    }  

    // initialize resampling indices 
//    getBootID(); done in R
    bootID = NULL;


    // Initialize GSL rnd environment variables
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    // an mt19937 generator with a seed of 0
    rnd = gsl_rng_alloc(T);
    if (mmRef->reprand!=TRUE){
       struct timeval tv;  // seed generation based on time
       gettimeofday(&tv, 0);
       unsigned long mySeed=tv.tv_sec + tv.tv_usec;
       gsl_rng_set(rnd, mySeed);  // reset seed
    }

//    printf("Anova test initialized.\n");

}

AnovaTest::~AnovaTest(){
}

void AnovaTest::releaseTest()
{
    free(multstat);
    free(Pmultstat);
    free(dfDiff);

    gsl_matrix_free(statj);
    gsl_matrix_free(Pstatj);
// The above commented out s.t. the results are passed out to R

    unsigned int i;
    for ( i=0; i<nModels; i++ ){
        gsl_matrix_free(Hats[i].mat);
        gsl_matrix_free(Hats[i].SS);
        gsl_matrix_free(Hats[i].R);
	gsl_matrix_free(Hats[i].Res);
        gsl_matrix_free(Hats[i].Coef);
	gsl_matrix_free(Hats[i].X);
	gsl_matrix_free(Hats[i].Y);
	gsl_vector_free(Hats[i].sd);
    }
    gsl_vector_free(bStatj);

    if ( bootID != NULL )
       gsl_matrix_free(bootID);

    for (i=0; i<nModels-1; i++)
       gsl_permutation_free(sortid[i]);
    free(sortid);

    gsl_rng_free(rnd);

//    printf("Anova test released.\n");

}

/*
void AnovaTest::display(void)
{
    unsigned int i, j, k;
    printf("Anova Table (resampling under H0):\n");
    printf("Hypo\t Alter\t df\t TestStat\t P-value\n");
    for ( i=0; i<nModels-1; i++ ){
       printf("Model%u\t Model%u\t %u\t %.3f\t\t %.3f\n", (unsigned int)i+1, (unsigned int)i, dfDiff[i], multstat[i], Pmultstat[i]);
    }

    if (mmRef->punit!=NONE){
       printf("Univariate Tests:\n");
       for (k=0; k<3; k++){
           printf("Response:\t");    
           for (j=k*4; j<(k+1)*4; j++)
               printf("variable%u\t", (unsigned int)j);    	
           printf("\n");
           for ( i=0; i<nModels-1; i++ ){       
               printf("Model %u:\t", (unsigned int)i+1);
               for (j=k*4; j<(k+1)*4; j++)
                   printf("%.3f(%.3f)\t", gsl_matrix_get(statj, i, j), gsl_matrix_get(Pstatj, i, j));
               printf("\n");
            }
	   printf("\n");
        }
     }	
}
*/

int AnovaTest::resampTest(void)
{
//    printf("Start resampling test ...\n");
    unsigned int i, j, p, id;
    unsigned int maxiter=mmRef->nboot; 
    double hii, score;

    gsl_matrix *bX, *bY;
    bY = gsl_matrix_alloc(nRows, nVars);
    bX = gsl_matrix_alloc(nRows, nParam);

    // initialize permid
    unsigned int *permid=NULL;
    if ( bootID == NULL ) {
       if ( mmRef->resamp == PERMUTE ){
          permid = (unsigned int *)malloc(nRows*sizeof(unsigned int));
          for (i=0; i<nRows; i++)
              permid[i] = i;
    } }
    // resampling options 
    if (mmRef->resamp == CASEBOOT) {
       nSamp = 0;
       for (i=0; i<maxiter; i++) {
           for ( j=0; j<nRows; j++ ){
	       // resampling index
 	       if (bootID == NULL) 
	          id = gsl_rng_uniform_int(rnd, nRows);
               else 
	          id = (unsigned int) gsl_matrix_get(bootID, i, j);
               // resample Y and X
               gsl_vector_view Yj=gsl_matrix_row(Yref, id);
               gsl_matrix_set_row (bY, j, &Yj.vector);
               gsl_vector_view Xj=gsl_matrix_row(Xref, id);
               gsl_matrix_set_row (bX, j, &Xj.vector); 
	    }
           anovacase(bY, bX);
           nSamp++;
        }
    } 
    else if (mmRef->resamp == RESIBOOT) {
        nSamp = 0;
        for (i=0; i<maxiter; i++) {
          for (p=1; p<nModels; p++) { 
            if (mmRef->reprand!=TRUE) {
                GetRNGstate();
                printf("reprand==FALSE\n");
            }
            for (j=0; j<nRows; j++){
               // resampling index
 	       if (bootID == NULL) 
	          id = gsl_rng_uniform_int(rnd, nRows);
               else 
	          id = (unsigned int) gsl_matrix_get(bootID, i, j);
               // bootr by resampling resi=(Y-fit)
               gsl_vector_view Yj=gsl_matrix_row(Yref, id);
               gsl_vector_view Fj=gsl_matrix_row(Hats[p].Y, id);
               gsl_matrix_set_row (bY, j, &Yj.vector);
               gsl_vector_view bootr=gsl_matrix_row(bY, j);
               gsl_vector_sub (&bootr.vector, &Fj.vector);  
               if (mmRef->student==TRUE) {
                  hii = gsl_matrix_get(Hats[p].mat, id, id);
                  gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
               } 
               // bY = Y + bootr
               Yj=gsl_matrix_row(Hats[p].Y, j);
               gsl_vector_add (&bootr.vector, &Yj.vector);
	    } 
            if (mmRef->reprand!=TRUE) PutRNGstate();
            anovaresi(bY, p);
         }
        nSamp++;
    } }
   else if (mmRef->resamp == SCOREBOOT) {
       nSamp = 0;
       for (i=0; i<maxiter; i++) {
         for (p=1; p<nModels; p++) {
           for ( j=0; j<nRows; j++ ) {
               // random score
	       if ( bootID == NULL )
	          score = gsl_ran_ugaussian (rnd); 
	       else
	          score = (double)gsl_matrix_get(bootID, i, j);
               // bootr = (Y - fit)*score 
               gsl_vector_view Yj=gsl_matrix_row(Yref, j);
               gsl_vector_view Fj=gsl_matrix_row(Hats[p].Y, j);
               gsl_matrix_set_row (bY, j, &Yj.vector);
               gsl_vector_view bootr=gsl_matrix_row(bY, j);
               gsl_vector_sub (&bootr.vector, &Fj.vector); 
               if (mmRef->student==TRUE) {
                  hii = gsl_matrix_get(Hats[p].mat, j, j);
                  gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
               }
                // bY = Y + bootr
               gsl_vector_scale (&bootr.vector, score);
               gsl_vector_add (&bootr.vector, &Fj.vector);
 	   } 
          anovaresi(bY, p);
        } 
        nSamp++;
   } }
   else if ( mmRef->resamp == PERMUTE ) { 
       gsl_matrix_add_constant (Pstatj, 1.0); 
       for (p=0; p<nModels-1; p++)
           Pmultstat[p]=1.0;       // include itself
        nSamp = 1;
        for (i=0; i<maxiter-1; i++) { //999
            for (p=1; p<nModels; p++){ 
                if (bootID == NULL ) 
                    gsl_ran_shuffle(rnd, permid, nRows, sizeof(unsigned int));
             // get bootr by permuting resi:Y-fit
                for (j=0; j<nRows; j++){
 	            if (bootID == NULL) 
	               id = permid[j];
                    else 
	               id = (unsigned int) gsl_matrix_get(bootID, i, j);
                   // bootr by resampling resi=(Y-fit)
                    gsl_vector_view Yj=gsl_matrix_row(Yref, id);
                    gsl_vector_view Fj=gsl_matrix_row(Hats[p].Y, id);
                    gsl_matrix_set_row (bY, j, &Yj.vector);
                    gsl_vector_view bootr=gsl_matrix_row(bY, j);
                    gsl_vector_sub (&bootr.vector, &Fj.vector); 
                    if (mmRef->student==TRUE) {
                        hii = gsl_matrix_get(Hats[p].mat, id, id);
                        gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
                    }
                    // bY = Y + bootr
                    Yj=gsl_matrix_row(Hats[p].Y, j);
                    gsl_vector_add (&bootr.vector, &Yj.vector);
                 }
                 anovaresi(bY, p);
           }
           nSamp++;
       }      
   }
   else 
       GSL_ERROR("Invalid resampling option", GSL_EINVAL);

   // p-values 
   unsigned int sid, sid0;
   double *pj;  
   for (i=0; i<nModels-1; i++) { 
        Pmultstat[i]=(double) (Pmultstat[i]+1)/(nSamp+1); // adjusted with +1
        pj = gsl_matrix_ptr (Pstatj, i, 0);
        if ( mmRef->punit == FREESTEP ){ 
           for (j=1; j<nVars; j++){
               sid = gsl_permutation_get(sortid[i], j);
	       sid0 = gsl_permutation_get(sortid[i], j-1);
	       *(pj+sid)=MAX(*(pj+sid), *(pj+sid0)); 
	   }  
        }
        if ( mmRef->punit == STEPUP ){ 
           for (j=2; j<nVars; j++){
               sid = gsl_permutation_get(sortid[i], nVars-j);
	       sid0 = gsl_permutation_get(sortid[i], nVars-j+1);
	       *(pj+sid) = MIN(*(pj+sid), *(pj+sid0)); 
	   }  
        }
        for (j=0; j<nVars; j++)
            *(pj+j) = (double)(*(pj+j)+1)/(nSamp+1);  // adjusted with +1 
    }

   // free memory
   gsl_matrix_free(bX);
   gsl_matrix_free(bY);
   if (permid!=NULL) free(permid);

   return 0;

}

int AnovaTest::anovacase(gsl_matrix *bY, gsl_matrix *bX)
{
   unsigned int j;
   // if Y col is all zeros
   for ( j=0; j<nVars; j++ ){
       gsl_vector_view colj = gsl_matrix_column(bY, j);
       if ( gsl_vector_isnull(&colj.vector) == TRUE ) return GSL_ERANGE;
   }

   unsigned int i, hid, aid;
   double *sj, *pj, *bj;
   gsl_matrix *Z = gsl_matrix_alloc(nRows, nVars);
   gsl_matrix_memcpy(Z, bY);
   // Hats.X 
   for (i=0; i<nModels-1; i++){
      hid = i+1; aid = i;  
      gsl_vector_view ref1 = gsl_matrix_row(inRef, aid);
      subX(bX, &ref1.vector, Hats[aid].X);
      gsl_vector_view ref0 = gsl_matrix_row(inRef, hid);
      subX(bX, &ref0.vector, Hats[hid].X);
      //Y = X*coef
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,Hats[aid].X,Hats[aid].Coef,0.0,Z); 
      //Z = bY - Yhat;
      gsl_matrix_add (Z, bY);
      // calc teststats
      calcSS(Z, &(Hats[hid]), mmRef);
      calcSS(Z, &(Hats[aid]), mmRef);
      testStatCalc(&(Hats[hid]), &(Hats[aid]), mmRef, TRUE, &(bMultStat), bStatj);

      if (bMultStat >= multstat[i]) Pmultstat[i]++;
      sj = gsl_matrix_ptr (statj, i, 0);
      pj = gsl_matrix_ptr (Pstatj, i, 0);
      bj = gsl_vector_ptr (bStatj, 0);          
      calcAdjustP(mmRef->punit, nVars, bj, sj, pj, sortid[i]);
   }

  gsl_matrix_free(Z);

  return 0;
}

int AnovaTest::anovaresi(gsl_matrix *bY, const unsigned int i)
{
    unsigned int hid=i, aid = i-1;

    // count the right-hand tails
    calcSS(bY, &(Hats[aid]), mmRef);
    calcSS(bY, &(Hats[hid]), mmRef);
    testStatCalc(&(Hats[hid]), &(Hats[aid]), mmRef, TRUE, &(bMultStat), bStatj);

    // count data related to P-values
    if (bMultStat >= multstat[aid]) Pmultstat[aid]++;
    // get result ptr corresponds to model i
    double *sj = gsl_matrix_ptr (statj, aid, 0);
    double *pj = gsl_matrix_ptr (Pstatj, aid, 0);
    double *bj = gsl_vector_ptr (bStatj, 0);
    calcAdjustP(mmRef->punit, nVars, bj, sj, pj, sortid[aid]);    
       
   return 0;
}
