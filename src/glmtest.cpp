// hypothesis testing, including summary and anova
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// // 16-3-2012

#include "resampTest.h"
//#include "time.h"
#include "multi_thread.h"
// constructor
GlmTest::GlmTest(const mv_Method *tm) : tm(tm) {
  eps = tm->tol;
  //    eps = 1e-6;

  smryStat = NULL;
  Psmry = NULL;

  anovaStat = NULL;
  Panova = NULL;

  Xin = NULL;

  XBeta = NULL;
  Sigma = NULL;
  bootID = NULL;
  bootStore = NULL;

  // Prepared for geeCalc
  L = gsl_matrix_alloc(tm->nParam, tm->nParam);
  gsl_matrix_set_identity(L);
  //Rlambda = gsl_matrix_alloc(tm->nVars, tm->nVars);

  Wj = gsl_matrix_alloc(tm->nRows, tm->nRows);

  // Initialize GSL rnd environment variables
  const gsl_rng_type *T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  // an mt19937 generator with a seed of 0
  rnd = gsl_rng_alloc(T);
  if (tm->reprand != TRUE) {
    struct timeval tv; // seed generation based on time
    gettimeofday(&tv, 0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    gsl_rng_set(rnd, mySeed); // reset seed
  }

  if (tm->resamp == PERMUTE) {
    permid = (size_t *)malloc(tm->nRows * sizeof(size_t));
    for (size_t i = 0; i < tm->nRows; i++)
      permid[i] = i;
  } else
    permid = NULL;

  if (tm->resamp == MONTECARLO) {
    XBeta = gsl_matrix_alloc(tm->nRows, tm->nVars);
    Sigma = gsl_matrix_alloc(tm->nVars, tm->nVars);
  }
  aic = new double[tm->nVars];

	thread_num = 0;
}

GlmTest::~GlmTest(void) {}

void GlmTest::init_multthread(int num) {
	if (num <= 0) return;
	thread_num = num;
}

// cleanup
void GlmTest::releaseTest(void) {
  if (smryStat != NULL)
    gsl_matrix_free(smryStat);
  if (Psmry != NULL)
    gsl_matrix_free(Psmry);

  if (anovaStat != NULL)
    gsl_matrix_free(anovaStat);
  if (Panova != NULL)
    gsl_matrix_free(Panova);
  if (bootStore != NULL)
    gsl_matrix_free(bootStore);


  gsl_matrix_free(L);
  //gsl_matrix_free(Rlambda);
  gsl_matrix_free(Wj);

  gsl_rng_free(rnd);
  if (XBeta != NULL)
    gsl_matrix_free(XBeta);
  if (Sigma != NULL)
    gsl_matrix_free(Sigma);

  if (bootID != NULL)
    gsl_matrix_free(bootID);
  if (permid != NULL)
    free(permid);

  delete[] aic;
}

int GlmTest::summary(glm *fit) {

  gsl_matrix *Rlambda =  gsl_matrix_alloc(tm->nVars, tm->nVars);

  unsigned int k;
  unsigned int nRows = tm->nRows, nVars = tm->nVars, nParam = tm->nParam;
  unsigned int mtype = fit->mmRef->model - 1;
  // REDESIGN THIS
  // Declare null models
  PoissonGlm pNull(fit->mmRef), pAlt(fit->mmRef);
  BinGlm binNull(fit->mmRef), binAlt(fit->mmRef);
  NBinGlm nbNull(fit->mmRef), nbAlt(fit->mmRef);
  GammaGlm gammaNull(fit->mmRef), gammaAlt(fit->mmRef);

  glm *PtrNull[4] = {&pNull, &nbNull, &binNull, &gammaNull};
  glm *PtrAlt[4] = {&pAlt, &nbAlt, &binAlt, &gammaAlt};

  gsl_vector_view teststat, unitstat;
  gsl_matrix_view L1;
  // To estimate initial Beta from PtrNull->Beta
  //  gsl_vector *ref=gsl_vector_alloc(nParam);
  //  gsl_matrix *BetaO=gsl_matrix_alloc(nParam, nVars);

  smryStat = gsl_matrix_alloc((nParam + 1), nVars + 1);
  Psmry = gsl_matrix_alloc((nParam + 1), nVars + 1);
  gsl_matrix_set_zero(Psmry);
  // initialize the design matrix for all hypo tests
  GrpMat *GrpXs = (GrpMat *)malloc((nParam + 2) * sizeof(GrpMat));
  GrpXs[0].matrix = gsl_matrix_alloc(nRows, nParam);
  gsl_matrix_memcpy(GrpXs[0].matrix, fit->Xref); // the alt X
  GrpXs[1].matrix = gsl_matrix_alloc(nRows, 1);  // overall test
  gsl_matrix_set_all(GrpXs[1].matrix, 1.0);

  // significance tests
  for (k = 2; k < nParam + 2; k++) {
    // if nParam = 1, ie just intercepts, this breaks I think this might be part
    // of issue #54
    GrpXs[k].matrix = gsl_matrix_alloc(nRows, nParam - 1);
    subX2(fit->Xref, k - 2, GrpXs[k].matrix);
  }
  // printf("\nmm->model %d, tm->test %d, tm->punit = %d, tm->corr = %d\n",
  //     fit->mmRef->model, tm->test, tm->punit, tm->corr);
  PtrNull[mtype]->initialGlm(fit->Yref, fit->Xref, fit->Oref, fit->Beta);
  PtrAlt[mtype]->initialGlm(fit->Yref, fit->Xref, fit->Oref, fit->Beta);


  // REDESIGN - make this a switch maybe
  // Calc test statistics
  if (tm->test == WALD) {
    // the overall test compares to mean
    teststat = gsl_matrix_row(smryStat, 0);
    L1 = gsl_matrix_submatrix(L, 1, 0, nParam - 1, nParam);
    double lambda = gsl_vector_get(tm->smry_lambda, 0);
    GetR(fit->Res, tm->corr, lambda, Rlambda);
    // Errors here if we have less than 2 predictors
    GeeWald(fit, &L1.matrix, &teststat.vector, Rlambda);
    // the significance test
    for (k = 2; k < nParam + 2; k++) {
      teststat = gsl_matrix_row(smryStat, k - 1);
      L1 = gsl_matrix_submatrix(L, k - 2, 0, 1, nParam);
      GeeWald(fit, &L1.matrix, &teststat.vector, Rlambda);
    }
  } else if (tm->test == SCORE) {
    for (k = 1; k < nParam + 2; k++) {
      teststat = gsl_matrix_row(smryStat, k - 1);
      PtrNull[mtype]->regression(fit->Yref, GrpXs[k].matrix, fit->Oref, NULL);
      double lambda = gsl_vector_get(tm->smry_lambda, k);
      GetR(PtrNull[mtype]->Res, tm->corr, lambda, Rlambda);
      GeeScore(GrpXs[0].matrix, PtrNull[mtype], &teststat.vector, Rlambda);
    }
  } else {
    for (k = 1; k < nParam + 2; k++) {
      teststat = gsl_matrix_row(smryStat, k - 1);
      PtrNull[mtype]->regression(fit->Yref, GrpXs[k].matrix, fit->Oref, NULL);
      GeeLR(fit, PtrNull[mtype], &teststat.vector); // works better
    }
  }
  // sort id if the unitvaraite test is free step-down
  gsl_permutation **sortid;
  sortid = (gsl_permutation **)malloc((nParam + 1) * sizeof(gsl_permutation *));
  for (k = 0; k < (nParam + 1); k++) {
    teststat = gsl_matrix_row(smryStat, k);
    unitstat = gsl_vector_subvector(&teststat.vector, 1, nVars);
    sortid[k] = gsl_permutation_alloc(nVars);
    gsl_sort_vector_index(sortid[k], &unitstat.vector);
    gsl_permutation_reverse(sortid[k]); // rearrange in descending order
  }

  if (tm->resamp == MONTECARLO) {
    double lambda = gsl_vector_get(tm->smry_lambda, 0);
    GetR(fit->Res, tm->corr, lambda, Sigma);
    setMonteCarlo(fit, XBeta, Sigma);
  }

  nSamp = 0;
  double *suj, *buj, *puj;
  gsl_matrix *bStat = gsl_matrix_alloc((nParam + 1), nVars + 1);
  gsl_matrix_set_zero(bStat);
  gsl_matrix *bY = gsl_matrix_alloc(nRows, nVars);
  gsl_matrix *bO = gsl_matrix_alloc(nRows, nVars);
  gsl_matrix_memcpy(bO, fit->Eta);
  double diff, timelast = 0;
  clock_t clk_start = clock();
  for (unsigned int i = 0; i < tm->nboot; i++) {
    nSamp++;
    if (tm->resamp == CASEBOOT)
      resampSmryCase(fit, bY, GrpXs, bO, i);
    else
      resampNonCase(fit, bY, i);

    if (tm->test == WALD) {
      PtrAlt[mtype]->regression(bY, GrpXs[0].matrix, bO, NULL);
      // the overall test compares to mean
      teststat = gsl_matrix_row(bStat, 0);
      L1 = gsl_matrix_submatrix(L, 1, 0, nParam - 1, nParam);
      double lambda = gsl_vector_get(tm->smry_lambda, 0);
      GetR(PtrAlt[mtype]->Res, tm->corr, lambda, Rlambda);
      GeeWald(PtrAlt[mtype], &L1.matrix, &teststat.vector, Rlambda);
      // the significance test
      for (k = 2; k < nParam + 2; k++) {
        teststat = gsl_matrix_row(bStat, k - 1);
        L1 = gsl_matrix_submatrix(L, k - 2, 0, 1, nParam);
        GeeWald(PtrAlt[mtype], &L1.matrix, &teststat.vector, Rlambda);
      }
    } else if (tm->test == SCORE) {
      for (k = 1; k < nParam + 2; k++) {
        teststat = gsl_matrix_row(bStat, k - 1);
        PtrNull[mtype]->regression(bY, GrpXs[k].matrix, bO, NULL);
        double lambda = gsl_vector_get(tm->smry_lambda, k);
        GetR(PtrNull[mtype]->Res, tm->corr, lambda, Rlambda);
        GeeScore(GrpXs[0].matrix, PtrNull[mtype], &teststat.vector, Rlambda);
      }
    } else { // use single bAlt estimate works better
      PtrAlt[mtype]->regression(bY, GrpXs[0].matrix, bO, NULL);
      for (k = 1; k < nParam + 2; k++) {
        teststat = gsl_matrix_row(bStat, k - 1);
        PtrNull[mtype]->regression(bY, GrpXs[k].matrix, bO, NULL);
        GeeLR(PtrAlt[mtype], PtrNull[mtype], &teststat.vector);
      }
    }
    for (k = 0; k < (nParam + 1); k++) {
      buj = gsl_matrix_ptr(bStat, k, 0);
      suj = gsl_matrix_ptr(smryStat, k, 0);
      puj = gsl_matrix_ptr(Psmry, k, 0);
      if (*buj >= *suj)
        *puj = *puj + 1;
      calcAdjustP(tm->punit, nVars, buj + 1, suj + 1, puj + 1, sortid[k]);
    } // end for j loop
    // Prompts
    if ((tm->showtime == TRUE) & (i % 100 == 0)) {
      diff = (float)(clock() - clk_start) / (float)CLOCKS_PER_SEC;
      timelast += (double)diff / 60;
      printf("\tResampling run %d finished. Time elapsed: %.2f min ...\n", i,
             timelast);
      clk_start = clock();
    }
  } // end for i loop

  // ========= Get P-values ========= //
  if (tm->punit == FREESTEP) {
    for (k = 0; k < (nParam + 1); k++) {
      puj = gsl_matrix_ptr(Psmry, k, 1);
      reinforceP(puj, nVars, sortid[k]);
    }
  }
  // p = (#exceeding observed stat + 1)/(#nboot+1)
  //    printf("tm->nboot=%d, nSamp=%d\n", tm->nboot, nSamp);
  gsl_matrix_add_constant(Psmry, 1.0);
  gsl_matrix_scale(Psmry, (double)1.0 / (nSamp + 1));

  for (k = 0; k < nVars; k++)
    aic[k] = -fit->ll[k] + 2 * (nParam + 1);

  // === release memory ==== //
  PtrAlt[mtype]->releaseGlm();
  if (tm->test != WALD)
    PtrNull[mtype]->releaseGlm();
  gsl_matrix_free(bStat);
  gsl_matrix_free(bY);
  gsl_matrix_free(bO);

  for (k = 0; k < nParam + 1; k++)
    if (sortid[k] != NULL)
      gsl_permutation_free(sortid[k]);
  free(sortid);

  if (GrpXs != NULL) {
    for (unsigned int k = 0; k < nParam + 2; k++)
      if (GrpXs[k].matrix != NULL)
        gsl_matrix_free(GrpXs[k].matrix);
    free(GrpXs);
  }
  gsl_matrix_free(Rlambda);

  return SUCCESS;
}

int GlmTest::anova_mt(glm *fit, gsl_matrix *isXvarIn){
  // Assume the models have been already sorted (in R)
  gsl_matrix *Rlambda =  gsl_matrix_alloc(tm->nVars, tm->nVars);
	Xin = isXvarIn;
  nModels = Xin->size1;
  double *rdf = new double[nModels];
  unsigned int nP, i, j;
  unsigned int ID0, ID1, nP0, nP1;
  unsigned int nRows = tm->nRows, nVars = tm->nVars, nParam = tm->nParam;
  unsigned int mtype = fit->mmRef->model - 1;

  dfDiff = new unsigned int[nModels - 1];
  anovaStat = gsl_matrix_alloc((nModels - 1), nVars + 1);
  Panova = gsl_matrix_alloc((nModels - 1), nVars + 1);
  bootStore = gsl_matrix_alloc(tm->nboot, nVars + 1);
	
  gsl_matrix_set_zero(anovaStat);
  gsl_matrix_set_zero(Panova);
  PoissonGlm pNull(fit->mmRef), pAlt(fit->mmRef);
  BinGlm binNull(fit->mmRef), binAlt(fit->mmRef);
  NBinGlm nbNull(fit->mmRef), nbAlt(fit->mmRef);
  GammaGlm gammaNull(fit->mmRef), gammaAlt(fit->mmRef);
  glm *PtrNull[4] = {&pNull, &nbNull, &binNull, &gammaNull};
  glm *PtrAlt[4] = {&pAlt, &nbAlt, &binAlt, &gammaAlt};

  
	glm *bNull[thread_num][4];
  glm *bAlt[thread_num][4];

	for (int k = 0; k < thread_num; ++k) {
	  bNull[k][0] = new PoissonGlm(fit->mmRef);
	  bNull[k][1] = new NBinGlm(fit->mmRef);
	  bNull[k][2] = new BinGlm(fit->mmRef);
	  bNull[k][3] = new GammaGlm(fit->mmRef);
	  bAlt[k][0] = new PoissonGlm(fit->mmRef);
	  bAlt[k][1] = new NBinGlm(fit->mmRef);
	  bAlt[k][2] = new BinGlm(fit->mmRef);
	  bAlt[k][3] = new GammaGlm(fit->mmRef);
	}
  
  double *suj, *buj, *puj;
  gsl_vector_view teststat, unitstat, ref1, ref0;
  gsl_matrix** X0 =  (gsl_matrix**) malloc(thread_num* sizeof (gsl_matrix*));
  gsl_matrix** X1 =  (gsl_matrix**) malloc(thread_num* sizeof (gsl_matrix*));
	gsl_matrix *L1 = NULL, *tmp1 = NULL;
	gsl_matrix **BetaO = (gsl_matrix**) malloc(thread_num*sizeof(gsl_matrix*));
  gsl_matrix **bO =(gsl_matrix**) malloc(thread_num*sizeof(gsl_matrix*)); 
	gsl_matrix** bY =(gsl_matrix**) malloc(thread_num*sizeof(gsl_matrix*)); 
	for (int k = 0; k < thread_num; ++k) {	
	  bY[k] = gsl_matrix_alloc(nRows, nVars);
		bO[k] = gsl_matrix_alloc(nRows, nVars);
	}
	for (int k = 0; k < thread_num; ++k) {
		BetaO[k] = NULL;
		X0[k] = NULL;
		X1[k] = NULL;
	}
  gsl_permutation *sortid = NULL;
  if (tm->punit == FREESTEP)
    sortid = gsl_permutation_alloc(nVars);

  // ======= Fit the (first) Alt model =========//
  for (i = 0; i < nModels; i++) {
    nP = 0;
    for (int k = 0; k < nParam; k++)
      if (gsl_matrix_get(Xin, i, k) != FALSE)
        nP++;
    rdf[i] = nRows - nP;
  }
	//init glm
	PtrNull[mtype]->initialGlm(fit->Yref, fit->Xref, fit->Oref, fit->Beta);
  PtrAlt[mtype]->initialGlm(fit->Yref, fit->Xref, fit->Oref, fit->Beta);

  for (i = 1; i < nModels; i++) {
    // ======= Fit the Null model =========//
    ID0 = i;
    ID1 = i - 1;
    nP0 = nRows - (unsigned int)rdf[ID0];
    nP1 = nRows - (unsigned int)rdf[ID1];
    
		// Degrees of freedom
    dfDiff[i - 1] = nP1 - nP0;

    ref1 = gsl_matrix_row(Xin, ID1);
    ref0 = gsl_matrix_row(Xin, ID0);
		for (int k = 0; k < thread_num; ++k) {
    	X0[k] = gsl_matrix_alloc(nRows, nP0);
    	subX(fit->Xref, &ref0.vector, X0[k]);
    	X1[k] = gsl_matrix_alloc(nRows, nP1);
    	subX(fit->Xref, &ref1.vector, X1[k]);
		}
    // ======= Get multivariate test statistics =======//
    // Estimate shrinkage parametr only once under H1
    // See "FW: Doubts R package "mvabund" (12/14/11)
    teststat = gsl_matrix_row(anovaStat, (i - 1));
		PtrNull[mtype]->regression(NULL, X0[0], fit->Oref, NULL);
		if (tm->test == SCORE) {
			double lambda = gsl_vector_get(tm->anova_lambda, ID0);
			GetR(PtrNull[mtype]->Res, tm->corr, lambda, Rlambda);
			GeeScore(X1[0], PtrNull[mtype], &teststat.vector, Rlambda);
		} else if (tm->test == WALD) {
			PtrAlt[mtype]->regression(NULL, X1[0], fit->Oref, NULL);
			L1 = gsl_matrix_alloc(nP1 - nP0, nP1);
			tmp1 = gsl_matrix_alloc(nParam, nP1);
			subX(L, &ref1.vector, tmp1);
			subXrow1(tmp1, &ref0.vector, &ref1.vector, L1);
			double lambda = gsl_vector_get(tm->anova_lambda, ID1);
			GetR(PtrAlt[mtype]->Res, tm->corr, lambda, Rlambda);
			GeeWald(PtrAlt[mtype], L1, &teststat.vector, Rlambda);
		} else { // test is LR
			BetaO[0] = gsl_matrix_alloc(nP1, nVars);
			addXrow2(PtrNull[mtype]->Beta, &ref1.vector, BetaO[0]);
			PtrAlt[mtype]->regression(NULL, X1[0], fit->Oref, BetaO[0]);
			GeeLR(PtrAlt[mtype], PtrNull[mtype], &teststat.vector);
		}
		if (tm->resamp == MONTECARLO) {
			double lambda = gsl_vector_get(tm->anova_lambda, ID0);
			GetR(fit->Res, tm->corr, lambda, Sigma);
			setMonteCarlo(PtrNull[mtype], XBeta, Sigma);
		}

		if (BetaO[0] != NULL) {
			for (int k = 1; k < thread_num; ++k) {
				BetaO[k] = gsl_matrix_alloc(nP1, nVars);
				gsl_matrix_memcpy(BetaO[k], BetaO[0]);	
			}
		}
    // ======= Get univariate test statistics =======//
  	for (int k = 0; k < thread_num; ++k) {	
			bAlt[k][mtype]->initialGlm(fit->Yref, X1[k], fit->Oref, NULL);
  		bNull[k][mtype]->initialGlm(fit->Yref, X0[k], fit->Oref, NULL);
		}
		
		//int stacksize = 1024 * 128;
		//pthread_attr_t attr;
		//int attr_ret = pthread_attr_init(&attr); 
		//attr_ret = pthread_attr_setstacksize(&attr, stacksize);
	  nSamp = tm->nboot;

		pthread_t* threadid = new pthread_t[thread_num]; 
		int thread_task = tm->nboot / thread_num + 1; 
		int last = 0;
		anovaboot* thread_data[thread_num];
		for (int k = 0; k < thread_num; ++k) {
			threadid[k] = k;
			thread_data[k] = (anovaboot*) malloc( 1* sizeof (anovaboot)); 
			anovaboot* th_data = thread_data[k];
			th_data->fit = fit;
			th_data->PtrAlt = PtrAlt[mtype];
			th_data->PtrNull = PtrNull[mtype];
			th_data->ref0 = &ref0.vector;;
			th_data->ref1 = &ref1.vector;
			th_data->L1 = L1;
			th_data->ID0 = ID0;
			th_data->ID1 = ID1;
			th_data->tm = (mv_Method*) tm;
			th_data->bAlt = bAlt[k][mtype];
			th_data->bNull = bNull[k][mtype];
			th_data->bY = bY[k];
			th_data->X1 = X1[k];
			th_data->X0 = X0[k];
			th_data->bO = bO[k];
			th_data->BetaO = BetaO[k];
			th_data->bootStore = bootStore;
		  
			th_data->bootID = bootID;
			th_data->rnd = rnd;
			th_data->XBeta = XBeta;
			th_data->Sigma = Sigma;
				
			th_data->start_counter = last;
			if (last + thread_task < tm->nboot) {
				th_data->loop_cnt = thread_task;
			} else {
				th_data->loop_cnt = tm->nboot - last;	
			}
			last += th_data->loop_cnt;
    	pthread_create(&threadid[k], NULL, &anovaboot_mt, th_data);
			//pthread_join (threadid[k], NULL);
		}
		for (int k = 0; k < thread_num; ++k) {
			pthread_join (threadid[k], NULL);
		}
		delete threadid;
 		//pthread_attr_destroy(&attr);
		for (int k = 0; k < thread_num; ++k) free(thread_data[k]);
    
		if (tm->punit == FREESTEP) {
      unitstat = gsl_vector_subvector(&teststat.vector, 1, nVars);
      gsl_sort_vector_index(sortid, &unitstat.vector);
      gsl_permutation_reverse(sortid);
    }
		gsl_vector* buv = gsl_vector_alloc(nVars + 1);

		for (int k = 0; k < tm->nboot; ++k) {
			gsl_matrix_get_row(buv, bootStore, k);
      //buj = gsl_matrix_ptr(bootStore, k, 0);
			buj = gsl_vector_ptr(buv, 0);
      suj = gsl_matrix_ptr(anovaStat, i - 1, 0);
      puj = gsl_matrix_ptr(Panova, i - 1, 0);
      if (*(buj) > (*(suj)-1e-8))
        *puj = *puj + 1;
      // ------ get univariate counts ---------//
      calcAdjustP(tm->punit, nVars, buj + 1, suj + 1, puj + 1, sortid);
		}
		gsl_vector_free(buv);
    // ========= get p-values ======== //
    if (tm->punit == FREESTEP) {
      puj = gsl_matrix_ptr(Panova, i - 1, 1);
      reinforceP(puj, nVars, sortid);
    }

    //  end for i loop
		
		for (int k = 0; k < thread_num; ++k) {
    	if (BetaO[k] != NULL) gsl_matrix_free(BetaO[k]);
    	if (X0[k] != NULL) gsl_matrix_free(X0[k]);
    	if (X1[k] != NULL) gsl_matrix_free(X1[k]);
		}

    if (tm->test == WALD) {
      if (L1 != NULL) gsl_matrix_free(L1);
      if (tmp1 != NULL) gsl_matrix_free(tmp1);
    }
  } // end i for loop  and test for loop
  
  // p = (#exceeding observed stat + 1)/(#nboot+1)
  gsl_matrix_add_constant(Panova, 1.0);
  gsl_matrix_scale(Panova, (double)1 / (nSamp + 1.0));
  
	//free resources
	for (int i = 0; i < thread_num; ++i) {
		bAlt[i][mtype]->releaseGlm();
		delete bAlt[i][0];
	  delete bAlt[i][1];
	  delete bAlt[i][2];
	  delete bAlt[i][3];
	 	bNull[i][mtype]->releaseGlm(); 
		delete bNull[i][0];
	  delete bNull[i][1];
	  delete bNull[i][2];
	  delete bNull[i][3];
	}

	PtrAlt[mtype]->releaseGlm();	  
  if (tm->test != WALD) {
   	PtrNull[mtype]->releaseGlm();
  }
  delete[] rdf;
  if (sortid != NULL)
    gsl_permutation_free(sortid);
  for (int i = 0; i < thread_num; ++i) {
  	gsl_matrix_free(bY[i]);
		if (bO[i] != NULL)	gsl_matrix_free(bO[i]);
 	}
  free(X0);
	free(X1);
	free(bY);
	free(BetaO);
	free(bO);
  /*for (j = 0; j < nVars + 1; j++) {
    printf("[Response %d]:", (unsigned int)j);
    for (i = 0; i < nModels - 1; i++)
      printf("\t%.3f (%.3f)", gsl_matrix_get(anovaStat, i, j),
             gsl_matrix_get(Panova, i, j));
    printf("\n");
  }*/
  gsl_matrix_free(Rlambda);
	//DisplayMatrix mat = DisplayMatrix(bootStore, "bootstore");
  //mat.Print();
	return SUCCESS;
}

int GlmTest::anova(glm *fit, gsl_matrix *isXvarIn) {
  // Assume the models have been already sorted (in R)
  gsl_matrix *Rlambda =  gsl_matrix_alloc(tm->nVars, tm->nVars);

  Xin = isXvarIn;
  nModels = Xin->size1;
  double *rdf = new double[nModels];
  unsigned int nP, i, j, k;
  unsigned int ID0, ID1, nP0, nP1;
  unsigned int nRows = tm->nRows, nVars = tm->nVars, nParam = tm->nParam;
  unsigned int mtype = fit->mmRef->model - 1;

  dfDiff = new unsigned int[nModels - 1];
  anovaStat = gsl_matrix_alloc((nModels - 1), nVars + 1);
  Panova = gsl_matrix_alloc((nModels - 1), nVars + 1);
  gsl_vector *bStat = gsl_vector_alloc(nVars + 1);
  bootStore = gsl_matrix_alloc(tm->nboot, nVars + 1);

  gsl_matrix_set_zero(anovaStat);
  gsl_matrix_set_zero(Panova);
  gsl_vector_set_zero(bStat);

  // There has to be a better way to do this
  PoissonGlm pNull(fit->mmRef), pAlt(fit->mmRef);
  BinGlm binNull(fit->mmRef), binAlt(fit->mmRef);
  NBinGlm nbNull(fit->mmRef), nbAlt(fit->mmRef);
  GammaGlm gammaNull(fit->mmRef), gammaAlt(fit->mmRef);

  PoissonGlm pNullb(fit->mmRef), pAltb(fit->mmRef);
  BinGlm binNullb(fit->mmRef), binAltb(fit->mmRef);
  NBinGlm nbNullb(fit->mmRef), nbAltb(fit->mmRef);
  GammaGlm gammaNullb(fit->mmRef), gammaAltb(fit->mmRef);

  glm *PtrNull[4] = {&pNull, &nbNull, &binNull, &gammaNull};
  glm *PtrAlt[4] = {&pAlt, &nbAlt, &binAlt, &gammaAlt};
  glm *bNull[4] = {&pNullb, &nbNullb, &binNullb, &gammaNullb};
  glm *bAlt[4] = {&pAltb, &nbAltb, &binAltb, &gammaAltb};

  double *suj, *buj, *puj;
  gsl_vector_view teststat, unitstat, ref1, ref0;
  gsl_matrix *X0 = NULL, *X1 = NULL, *L1 = NULL, *tmp1 = NULL, *BetaO = NULL;
  gsl_matrix *bO = NULL, *bY = gsl_matrix_alloc(nRows, nVars);
  bO = gsl_matrix_alloc(nRows, nVars);

  gsl_permutation *sortid = NULL;
  if (tm->punit == FREESTEP)
    sortid = gsl_permutation_alloc(nVars);

  // ======= Fit the (first) Alt model =========//
  for (i = 0; i < nModels; i++) {
    nP = 0;
    for (k = 0; k < nParam; k++)
      if (gsl_matrix_get(Xin, i, k) != FALSE)
        nP++;
    rdf[i] = nRows - nP;
  }
	//init glm
  PtrNull[mtype]->initialGlm(fit->Yref, fit->Xref, fit->Oref, fit->Beta);
  PtrAlt[mtype]->initialGlm(fit->Yref, fit->Xref, fit->Oref, fit->Beta);

  for (i = 1; i < nModels; i++) {
    // ======= Fit the Null model =========//
    ID0 = i;
    ID1 = i - 1;
    nP0 = nRows - (unsigned int)rdf[ID0];
    nP1 = nRows - (unsigned int)rdf[ID1];
    
		// Degrees of freedom
    dfDiff[i - 1] = nP1 - nP0;

    ref1 = gsl_matrix_row(Xin, ID1);
    ref0 = gsl_matrix_row(Xin, ID0);
    X0 = gsl_matrix_alloc(nRows, nP0);
    subX(fit->Xref, &ref0.vector, X0);
    X1 = gsl_matrix_alloc(nRows, nP1);
    subX(fit->Xref, &ref1.vector, X1);

    // ======= Get multivariate test statistics =======//
    // Estimate shrinkage parametr only once under H1
    // See "FW: Doubts R package "mvabund" (12/14/11)
    teststat = gsl_matrix_row(anovaStat, (i - 1));
    PtrNull[mtype]->regression(NULL, X0, fit->Oref, NULL);
    if (tm->test == SCORE) {
      double lambda = gsl_vector_get(tm->anova_lambda, ID0);
      GetR(PtrNull[mtype]->Res, tm->corr, lambda, Rlambda);
      GeeScore(X1, PtrNull[mtype], &teststat.vector, Rlambda);
    } else if (tm->test == WALD) {
      PtrAlt[mtype]->regression(NULL, X1, fit->Oref, NULL);
      L1 = gsl_matrix_alloc(nP1 - nP0, nP1);
      tmp1 = gsl_matrix_alloc(nParam, nP1);
      subX(L, &ref1.vector, tmp1);
      subXrow1(tmp1, &ref0.vector, &ref1.vector, L1);
      double lambda = gsl_vector_get(tm->anova_lambda, ID1);
      GetR(PtrAlt[mtype]->Res, tm->corr, lambda, Rlambda);
      GeeWald(PtrAlt[mtype], L1, &teststat.vector, Rlambda);
    } else { // test is LR
      BetaO = gsl_matrix_alloc(nP1, nVars);
      addXrow2(PtrNull[mtype]->Beta, &ref1.vector, BetaO);
      PtrAlt[mtype]->regression(NULL, X1, fit->Oref, BetaO);
      GeeLR(PtrAlt[mtype], PtrNull[mtype], &teststat.vector);
    }

    if (tm->resamp == MONTECARLO) {
      double lambda = gsl_vector_get(tm->anova_lambda, ID0);
      GetR(fit->Res, tm->corr, lambda, Sigma);
      setMonteCarlo(PtrNull[mtype], XBeta, Sigma);
    }

  	
		bAlt[mtype]->initialGlm(fit->Yref, X1, fit->Oref, NULL);
  	bNull[mtype]->initialGlm(fit->Yref, X0, fit->Oref, NULL);

    // ======= Get resampling distribution under H0 ===== //
    double dif, timelast = 0;
    clock_t clk_start = clock();
    if (tm->showtime == TRUE)
      printf("Resampling begins for test %d.\n", i);
    for (j = 0; j < tm->nboot; j++) {
      //            printf("simu %d :", j);
      gsl_vector_set_zero(bStat);
      if (tm->resamp == CASEBOOT) {
        resampAnovaCase(PtrAlt[mtype], bY, X1, bO, j);
        subX(X1, &ref0.vector, X0);
      } else {
        resampNonCase(PtrNull[mtype], bY, j);
        gsl_matrix_memcpy(bO, fit->Oref);
      }

      if (tm->test == WALD) {
				if (tm->resamp == CASEBOOT) {
        	bAlt[mtype]->regression(bY, X1, bO, NULL);
				} else {
			   	bAlt[mtype]->regression(bY, NULL, bO, NULL);
				}
				double lambda = gsl_vector_get(tm->anova_lambda, ID1);
        GetR(bAlt[mtype]->Res, tm->corr, lambda, Rlambda);
        GeeWald(bAlt[mtype], L1, bStat, Rlambda);
      } else if (tm->test == SCORE) {
				if (tm->resamp == CASEBOOT) {
        	bNull[mtype]->regression(bY, X0, bO, NULL);
				} else {
			   	bNull[mtype]->regression(bY, NULL, bO, NULL);
				}
				double lambda = gsl_vector_get(tm->anova_lambda, ID0);
        GetR(bNull[mtype]->Res, tm->corr, lambda, Rlambda);
        GeeScore(X1, bNull[mtype], bStat, Rlambda);
      } else {
				if (tm->resamp == CASEBOOT) {
        	bNull[mtype]->regression(bY, X0, bO, NULL);
				} else {
			   	bNull[mtype]->regression(bY, NULL, bO, NULL);
				}
				addXrow2(bNull[mtype]->Beta, &ref1.vector, BetaO);
				if (tm->resamp == CASEBOOT) {
        	bAlt[mtype]->regression(bY, X1, bO, BetaO);
				} else {
			   	bAlt[mtype]->regression(bY, NULL, bO, BetaO);
				}
				GeeLR(bAlt[mtype], bNull[mtype], bStat);
      }
      gsl_matrix_set_row(bootStore, j, bStat);
      // ----- get multivariate counts ------- //
      if ((tm->showtime == TRUE) & (j % 100 == 0)) {
        dif = (float)(clock() - clk_start) / (float)CLOCKS_PER_SEC;
        timelast += (double)dif / 60;
        printf("\tResampling run %d finished. Time elapsed: %.2f minutes...\n",
               j, timelast);
        clk_start = clock();
      }
		}
    
		
		// ======= Get univariate test statistics =======//
    if (tm->punit == FREESTEP) {
      unitstat = gsl_vector_subvector(&teststat.vector, 1, nVars);
      gsl_sort_vector_index(sortid, &unitstat.vector);
      gsl_permutation_reverse(sortid);
    }

		//do statistics
		for (int k = 0; k < tm->nboot; k++){
			buj = gsl_matrix_ptr(bootStore, k, 0);
      suj = gsl_matrix_ptr(anovaStat, i - 1, 0);
      puj = gsl_matrix_ptr(Panova, i - 1, 0);
      if (*(buj) > (*(suj)-1e-8))
        *puj = *puj + 1;
      // ------ get univariate counts ---------//
      calcAdjustP(tm->punit, nVars, buj + 1, suj + 1, puj + 1, sortid);
      // Prompts
    } // end j for loop

    // ========= get p-values ======== //
    if (tm->punit == FREESTEP) {
      puj = gsl_matrix_ptr(Panova, i - 1, 1);
      reinforceP(puj, nVars, sortid);
    }
    //  end for i loop

    if (BetaO != NULL)
      gsl_matrix_free(BetaO);
    if (X0 != NULL)
      gsl_matrix_free(X0);
    if (X1 != NULL)
      gsl_matrix_free(X1);
    if (tm->test == WALD) {
      if (L1 != NULL)
        gsl_matrix_free(L1);
      if (tmp1 != NULL)
        gsl_matrix_free(tmp1);
    }
  } // end i for loop  and test for loop
  
	nSamp = tm->nboot;
  // p = (#exceeding observed stat + 1)/(#nboot+1)
  gsl_matrix_add_constant(Panova, 1.0);
  gsl_matrix_scale(Panova, (double)1 / (nSamp + 1.0));

  bAlt[mtype]->releaseGlm();
  PtrAlt[mtype]->releaseGlm();
  if (tm->test != WALD) {
    bNull[mtype]->releaseGlm();
    PtrNull[mtype]->releaseGlm();
  }
  delete[] rdf;
  if (sortid != NULL)
    gsl_permutation_free(sortid);
  gsl_vector_free(bStat);
  gsl_matrix_free(bY);
  if (bO != NULL)
    gsl_matrix_free(bO);
	gsl_matrix_free(Rlambda);
  return SUCCESS;
}

int GlmTest::GeeLR(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat) {
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

int GlmTest::GeeScore(gsl_matrix *X1, glm *PtrNull, gsl_vector *teststat, gsl_matrix* rlambda) {
  gsl_set_error_handler_off();

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
int GlmTest::GeeWald(glm *Alt, gsl_matrix *LL, gsl_vector *teststat, gsl_matrix* rlambda) {
  gsl_set_error_handler_off();

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

int GlmTest::resampSmryCase(glm *model, gsl_matrix *bT, GrpMat *GrpXs,
                            gsl_matrix *bO, unsigned int i) {
  gsl_set_error_handler_off();
  int status, isValid = TRUE;

  unsigned int j, k, id;
  gsl_vector_view yj, oj, xj;
  unsigned int nRows = tm->nRows, nParam = tm->nParam;
  gsl_matrix *tXX = gsl_matrix_alloc(nParam, nParam);

  while (isValid == TRUE) {
    // if all isSingular==TRUE
    for (j = 0; j < nRows; j++) {
      // resample Y, X, offsets accordingly
      if (bootID != NULL) {
        id = (unsigned int)gsl_matrix_get(bootID, i, j);
      } else {
        if (tm->reprand == TRUE) {
          id = (unsigned int)gsl_rng_uniform_int(rnd, nRows);
        } else {
          id = (unsigned int)nRows * Rf_runif(0, 1);
        }
      }
      xj = gsl_matrix_row(model->Xref, id);
      gsl_matrix_set_row(GrpXs[0].matrix, j, &xj.vector);
      yj = gsl_matrix_row(model->Yref, id);
      gsl_matrix_set_row(bT, j, &yj.vector);
      oj = gsl_matrix_row(model->Eta, id);
      gsl_matrix_set_row(bO, j, &oj.vector);
    }
    gsl_matrix_set_identity(tXX);
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, GrpXs[0].matrix, 0.0, tXX);
    status = gsl_linalg_cholesky_decomp(tXX);
    if (status != GSL_EDOM)
      break;
  }

  for (k = 2; k < nParam + 2; k++) {
    subX2(GrpXs[0].matrix, k - 2, GrpXs[k].matrix);
  }

  gsl_matrix_free(tXX);

  return SUCCESS;
}

int GlmTest::resampAnovaCase(glm *model, gsl_matrix *bT, gsl_matrix *bX,
                             gsl_matrix *bO, unsigned int i) {
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

int GlmTest::resampNonCase(glm *model, gsl_matrix *bT, unsigned int i) {
  unsigned int j, k, id;
  double bt, score, yij, mij;
  gsl_vector_view yj;
  unsigned int nRows = tm->nRows, nVars = tm->nVars;
  // to store Rf_unif
  //   gsl_vector *tmp = gsl_vector_alloc(nRows);
  //   gsl_permutation *vperm = gsl_permutation_alloc(nRows);
  double *tmp = (double *)malloc(nRows * sizeof(double));

  // note that residuals have got means subtracted
  switch (tm->resamp) {
  case RESIBOOT:
    for (j = 0; j < nRows; j++) {
      if (bootID != NULL)
        id = (unsigned int)gsl_matrix_get(bootID, i, j);
      else if (tm->reprand == TRUE)
        id = (unsigned int)gsl_rng_uniform_int(rnd, nRows);
      else
        id = (unsigned int)nRows * Rf_runif(0, 1);
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
    for (j = 0; j < nRows; j++) {
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

  return SUCCESS;
}

void GlmTest::displaySmry(glm *fit) {
  unsigned int i, j, k, nk;
  unsigned int nVars = tm->nVars, nParam = tm->nParam;
  const char *testname[3]              // string array, only pointers stored
      = {"sqrt(WALD)", "SCORE", "LR"}; // 2, 3, 4
                                       /*
                                           printf("\nSummary of fitting (resampling under H1):\n");
                                           printf("\n - Regression performance: \n");
                                           printf("\t\t phi\t AIC\t log-like\t");
                                           for (j=0; j<nVars; j++){
                                               printf("\n[Response %d]:\t", (unsigned int)j+1);
                                               printf("%.4f\t ", fit->phi[j]);
                                               printf("%.2f\t ", aic[j]);
                                               printf("%.2f", fit->ll[j]);
                                           }
                                       */
  if (tm->corr == SHRINK)
    displayvector(tm->smry_lambda, "\n Est. shrink.param in summary\n");

  // significance test
  nk = 1;
  k = tm->test - 2;
  printf("\n - Significance test (Pr>=%s):\n", testname[k]);
  if (tm->punit == FREESTEP)
    printf("\t (FREESTEP adjusted)\n");
  while (nk < nParam + 1) {
    printf("\t");
    for (i = nk; i < MIN(nk + 4, nParam + 1); i++)
      printf("\t [Explain %d] ", (unsigned int)i);
    printf("\n\t ");
    for (i = nk; i < MIN(nk + 4, nParam + 1); i++)
      printf(" %.3f (%.3f) \t", gsl_matrix_get(smryStat, i, 0),
             gsl_matrix_get(Psmry, i, 0));
    printf("\n\n");
    // Significance univariate tests
    if (tm->punit > NONE) {
      for (j = 1; j < nVars + 1; j++) {
        printf("[Response %d]:\t", (int)j);
        for (i = nk; i < MIN(nk + 4, nParam + 1); i++)
          printf("%.3f (%.3f)\t", gsl_matrix_get(smryStat, i, j),
                 gsl_matrix_get(Psmry, i, j));
        printf("\n");
      }
    }
    nk = i;
    printf("\n");
  }
  // Overall statistics
  printf("\n - Multivariate test (Pr>=%s): %.3f (%.3f)", testname[k],
         gsl_matrix_get(smryStat, 0, 0), gsl_matrix_get(Psmry, 0, 0));
  if (tm->punit == FREESTEP) {
    printf("\t (FREESTEP adjusted)\n");
    for (j = 1; j < nVars + 1; j++)
      printf("[Response %d]:\t%.3f (%.3f)\n", (int)j,
             gsl_matrix_get(smryStat, 0, j), gsl_matrix_get(Psmry, 0, j));
  }
  printf("\n ========================= \n");
}

void GlmTest::displayAnova(void) {
  unsigned int nVars = tm->nVars;
  unsigned int i, j;
  const char *testname[3]              // string array, only pointers stored
      = {"sqrt(WALD)", "SCORE", "LR"}; // 2, 3, 4

  displaymatrix(bootID, "bootID");

  printf("\n ========================= \n");
  printf("\nAnova Table (resampling under ");
  if (tm->resamp == CASEBOOT)
    printf("H1):\n");
  else
    printf("H0):\n");

  if (tm->corr == SHRINK)
    displayvector(tm->anova_lambda, "Est. shrink.param in anova");

  unsigned int test = tm->test - 2;
  printf("Hypo\t Alter\t dff\t %s\t  P-value \n", testname[test]);
  for (i = 0; i < nModels - 1; i++)
    printf("M%d\t M%d\t %d\t %.3f   %.3f\t\t \n", (int)i + 1, (int)i, dfDiff[i],
           gsl_matrix_get(anovaStat, i, 0), gsl_matrix_get(Panova, i, 0));

  if (tm->punit != NONE) {
    if (tm->punit == FREESTEP)
      printf("\nUnivariate Tests (FREESTEP adjusted):\n\t\t");
    else
      printf("\nUnivariate Tests:\n\t\t");
    for (i = 0; i < nModels - 1; i++)
      printf("\tM%d v. M%d\t", (unsigned int)i + 1, (unsigned int)i);
    printf("\n");

    for (j = 1; j < nVars + 1; j++) {
      printf("[Response %d]:", (unsigned int)j);
      for (i = 0; i < nModels - 1; i++)
        printf("\t%.3f (%.3f)", gsl_matrix_get(anovaStat, i, j),
               gsl_matrix_get(Panova, i, j));
      printf("\n");
    }
    printf("\n");
  }
}
