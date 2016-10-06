// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Interface between R and anova.cpp (Rcpp API >= 0.7.11)
//
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// Last modified: 20-April-2010
//
// Rcpp/RcppGSL changes by Dirk Eddelbuettel, July - August 2015


#include <RcppGSL.h>
#include "mvabund_types.h"
extern "C"{
#include "resampTest.h"
#include "time.h"
}

using namespace Rcpp;

// declare a dependency on the headers in the RcppGSL package;
// also activates plugin
// [[Rcpp::depends(RcppGSL)]]

// declare the function to be 'exported' to R
// [[Rcpp::export]]
List RtoAnovaCpp(const List & rparam,
                 RcppGSL::Matrix & Y,
                 RcppGSL::Matrix & X,
                 RcppGSL::Matrix & isXvarIn,
                 Rcpp::Nullable<RcppGSL::Matrix> & bID)
{
    // pass parameters
    mv_Method mm;	
//    mm.tol = as<double>(rparam["tol"]);
    mm.nboot = as<unsigned int>(rparam["nboot"]);
    mm.corr = as<unsigned int>(rparam["cor_type"]);
    mm.shrink_param = as<double>(rparam["shrink_param"]);
    mm.test = as<unsigned int>(rparam["test_type"]);
    mm.resamp = as<unsigned int>(rparam["resamp"]);
    mm.reprand = as<unsigned int>(rparam["reprand"]);
    mm.student = as<unsigned int>(rparam["studentize"]);
    mm.punit = as<unsigned int>(rparam["punit"]);
    mm.rsquare = as<unsigned int>(rparam["rsquare"]);

    unsigned int nRows = Y.nrow();
    unsigned int nModels = isXvarIn.nrow();

// initialize anova class
    AnovaTest anova(&mm, Y, X, isXvarIn);
	
// Resampling indices
    if (bID.isNotNull()) {
        RcppGSL::Matrix m(bID);
        mm.nboot = m.nrow();	   
        anova.bootID = m;
    }
//    else Rprintf("bID is null");

    // resampling test
    anova.resampTest();
//    anova.display();

    // Wrap the gsl objects with Rcpp, then Rcpp -> R
    List rs = List::create(
        Named("multstat" )= NumericVector(anova.multstat, anova.multstat+nModels-1),
        Named("Pmultstat")= NumericVector(anova.Pmultstat, anova.Pmultstat+nModels-1),
        Named("dfDiff"   )= NumericVector(anova.dfDiff, anova.dfDiff+nModels-1),
        Named("statj"    )= RcppGSL::Matrix(anova.statj),
        Named("Pstatj"   )= RcppGSL::Matrix(anova.Pstatj),
        Named("nSamp"    )= anova.nSamp
    );

    // clear objects
    anova.releaseTest(); 

    return rs;
}


// declare the function to be 'exported' to R
// [[Rcpp::export]]
List RtoGlmAnova(const List & sparam,
                 const List & rparam,
                 RcppGSL::Matrix & Y,
                 RcppGSL::Matrix & X,
                 RcppGSL::Matrix & O,
                 RcppGSL::Matrix & isXvarIn,
                 Rcpp::Nullable<RcppGSL::Matrix> & bID,
                 RcppGSL::Vector & lambda)
{
    // pass regression parameters
    reg_Method mm;
    mm.tol = as<double>(sparam["tol"]);
    mm.model = as<unsigned int>(sparam["regression"]);
    mm.speclink = as<unsigned int>(sparam["link"]);
    mm.estiMethod = as<unsigned int>(sparam["estimation"]);
    mm.varStab = as<unsigned int>(sparam["stablizer"]);   
    mm.n = as<unsigned int>(sparam["n"]);   
    mm.maxiter = as<unsigned int>(sparam["maxiter"]);   
    mm.maxiter2 = as<unsigned int>(sparam["maxiter2"]);   
    mm.warning = as<unsigned int>(sparam["warning"]);   

    // pass test parameters
    mv_Method tm;	
    tm.nboot = as<unsigned int>(rparam["nboot"]);
    tm.corr = as<unsigned int>(rparam["cor_type"]);
    tm.test = as<unsigned int>(rparam["test_type"]);
    tm.resamp = as<unsigned int>(rparam["resamp"]);
    tm.reprand = as<unsigned int>(rparam["reprand"]);
    tm.punit = as<unsigned int>(rparam["punit"]);
    tm.showtime = as<unsigned int>(rparam["showtime"]);
    tm.warning = as<unsigned int>(rparam["warning"]);   

    unsigned int nRows = Y.nrow();
    unsigned int nVars = Y.ncol();
    unsigned int nParam = X.ncol();
    unsigned int nModels = isXvarIn.nrow();
    unsigned int nLambda = lambda.size();
  
    tm.anova_lambda = gsl_vector_alloc(nLambda);
    gsl_vector_memcpy(tm.anova_lambda, lambda);

    tm.nRows = nRows;
    tm.nVars = nVars;
    tm.nParam = nParam;

    // do stuff	
    clock_t clk_start, clk_end;
    clk_start = clock();

    // glmfit
    PoissonGlm pfit(&mm);
    NBinGlm nbfit(&mm);
    BinGlm binfit(&mm);
    glm *glmPtr[3] = { &pfit, &nbfit, &binfit };
    unsigned int mtype = mm.model-1;
    glmPtr[mtype]->regression(Y, X, O, NULL);
//    glmPtr[mtype]->display();

    GlmTest myTest(&tm);
    // Resampling indices
    if ( bID.isNotNull() ) {
        RcppGSL::Matrix m(bID);
        tm.nboot = m.nrow();
        myTest.bootID = m;
    }  

    // resampling test
    myTest.anova(glmPtr[mtype],isXvarIn);
//    myTest.displayAnova();

    // Timing
    clk_end = clock();
    unsigned long int dif = floor((double)(clk_end - clk_start)/(double)(CLOCKS_PER_SEC));  
    unsigned int hours = floor((double)(dif/(double)3600));
    unsigned int min = floor((double)(dif%3600)/(double)60);
    unsigned int sec = dif%60;   
    if (tm.showtime>=TRUE)
        Rprintf("Time elapsed: %d hr %d min %d sec\n", hours, min, sec);

    // Wrap the gsl objects with Rcpp, then Rcpp -> R 
    RcppGSL::VectorView mul= gsl_matrix_column(myTest.anovaStat,0);
    RcppGSL::VectorView Pmul=gsl_matrix_column(myTest.Panova, 0);
    RcppGSL::MatrixView stat=gsl_matrix_submatrix(myTest.anovaStat,0,1,nModels-1,nVars);
    RcppGSL::MatrixView Pstat=gsl_matrix_submatrix(myTest.Panova,0,1,nModels-1,nVars);

    List rs = List::create(
         Named("multstat" )= mul,
         Named("Pmultstat")= Pmul,
         Named("statj"     )= stat,
         Named("Pstatj"    )= Pstat,
         Named("dfDiff"   )= NumericVector(myTest.dfDiff, myTest.dfDiff+nModels-1),
         Named("nSamp"    )= myTest.nSamp
    );

    // clear objects
    myTest.releaseTest();
    glmPtr[mtype]->releaseGlm();
    gsl_vector_free(tm.anova_lambda);

    return rs;
}


// declare the function to be 'exported' to R
// [[Rcpp::export]]
List RtoGlm(const List & rparam,
            RcppGSL::Matrix & Y,
            RcppGSL::Matrix & X,
            RcppGSL::Matrix & O)
{

    // pass regression parameters
    reg_Method mm;	
    mm.tol = as<double>(rparam["tol"]);
    mm.model = as<unsigned int>(rparam["regression"]);
    mm.speclink = as<unsigned int>(rparam["link"]);
    mm.estiMethod = as<unsigned int>(rparam["estimation"]);
    mm.varStab = as<unsigned int>(rparam["stablizer"]);
    mm.n = as<unsigned int>(rparam["n"]);
    mm.maxiter = as<unsigned int>(rparam["maxiter"]);
    mm.maxiter2 = as<unsigned int>(rparam["maxiter2"]);
    mm.warning = as<unsigned int>(rparam["warning"]);

    // do stuff	
    PoissonGlm pfit(&mm);
    BinGlm lfit(&mm);
    NBinGlm nbfit(&mm);
    glm *glmPtr[3] = { &pfit, &nbfit, &lfit};
    unsigned int mtype = mm.model-1;
    glmPtr[mtype]->regression(Y, X, O, NULL);
//    glmPtr[mtype]->display();
	
    // Wrap the glm object with Rcpp, then Rcpp -> R  
    List rs = List::create(
       Named("coefficients") = RcppGSL::Matrix(glmPtr[mtype]->Beta),
       Named("var.coefficients") = RcppGSL::Matrix(glmPtr[mtype]->varBeta),
       Named("fitted.values") = RcppGSL::Matrix(glmPtr[mtype]->Mu),
       Named("linear.predictor") = RcppGSL::Matrix(glmPtr[mtype]->Eta),
       Named("residuals") = RcppGSL::Matrix(glmPtr[mtype]->Res),
       Named("PIT.residuals") = RcppGSL::Matrix(glmPtr[mtype]->PitRes),
       Named("sqrt.1_Hii") = RcppGSL::Matrix(glmPtr[mtype]->sqrt1_Hii),
       Named("var.estimator") = RcppGSL::Matrix(glmPtr[mtype]->Var),
       Named("sqrt.weight") = RcppGSL::Matrix(glmPtr[mtype]->wHalf),
       Named("theta")=NumericVector(glmPtr[mtype]->theta, glmPtr[mtype]->theta+Y.ncol()),
       Named("two.loglike")=NumericVector(glmPtr[mtype]->ll, glmPtr[mtype]->ll+Y.ncol()),
       Named("aic")=NumericVector(glmPtr[mtype]->aic, glmPtr[mtype]->aic+Y.ncol()),
       Named("deviance")=NumericVector(glmPtr[mtype]->dev, glmPtr[mtype]->dev+Y.ncol()),
       Named("iter")=NumericVector(glmPtr[mtype]->iterconv, glmPtr[mtype]->iterconv+Y.ncol())
    );

    // clear objects
    glmPtr[mtype]->releaseGlm();

    return rs;
}

// declare the function to be 'exported' to R
// [[Rcpp::export]]
List RtoGlmSmry(const List & sparam,
                const List & rparam,
                RcppGSL::Matrix & Y,
                RcppGSL::Matrix & X,
                RcppGSL::Matrix & O,
                Rcpp::Nullable<RcppGSL::Matrix> & bID,
                RcppGSL::Vector & lambda)
{
    // Pass regression parameters
    reg_Method mm;
    mm.tol = as<double>(sparam["tol"]);
    mm.model = as<unsigned int>(sparam["regression"]);
    mm.speclink = as<unsigned int>(sparam["link"]);
    mm.estiMethod = as<unsigned int>(sparam["estimation"]);
    mm.varStab = as<unsigned int>(sparam["stablizer"]);
    mm.n = as<unsigned int>(sparam["n"]);
    mm.maxiter = as<unsigned int>(sparam["maxiter"]);
    mm.maxiter2 = as<unsigned int>(sparam["maxiter2"]);
    mm.warning = as<unsigned int>(sparam["warning"]);

    // pass test parameters
    mv_Method tm;
    tm.corr = as<unsigned int>(rparam["cor_type"]);
    tm.test = as<unsigned int>(rparam["test_type"]);
    tm.resamp = as<unsigned int>(rparam["resamp"]);
    tm.reprand = as<unsigned int>(rparam["reprand"]);
    tm.punit = as<unsigned int>(rparam["punit"]);
    tm.nboot = as<unsigned int>(rparam["nboot"]);
    tm.showtime = as<unsigned int>(rparam["showtime"]);
    tm.warning = as<unsigned int>(rparam["warning"]);

    unsigned int nRows = Y.nrow();
    unsigned int nVars = Y.ncol();
    unsigned int nParam = X.ncol();
    unsigned int nLambda = lambda.size();

    tm.smry_lambda = gsl_vector_alloc(nLambda);
    gsl_vector_memcpy(tm.smry_lambda, lambda);

    tm.nRows = nRows;
    tm.nVars = nVars;
    tm.nParam = nParam;

    // do stuff	
    clock_t clk_start, clk_end;
    clk_start = clock();

    // Glm fit
    PoissonGlm pfit(&mm);
    BinGlm lfit(&mm);
    NBinGlm nbfit(&mm);
    glm *glmPtr[3] = { &pfit, &nbfit, &lfit };
    unsigned int mtype = mm.model-1;
    glmPtr[mtype]->regression(Y, X, O, NULL);
//    glmPtr[mtype]->display();

    GlmTest myTest(&tm);    
    // Resampling indices
    if ( bID.isNotNull() ) {
        RcppGSL::Matrix m(bID);
        tm.nboot = m.nrow();        
        myTest.bootID = m;
    }
    // resampling test
    myTest.summary(glmPtr[mtype]);
//    myTest.displaySmry(glmPtr[mtype]);

    // Timing
    clk_end = clock();
    unsigned long int dif = floor((double)(clk_end - clk_start)/(double)(CLOCKS_PER_SEC));
    unsigned int hours = floor((double)dif/(double)3600);
    unsigned int min = floor((double)(dif-hours*3600)/(double)60);
    unsigned int sec = dif - hours*3600 - min*60;
    if (tm.showtime>=TRUE)
       Rprintf("Time elapsed: %d hr %d min %d sec\n", hours, min, sec);

    // Wrap gsl vectors with Rcpp, then Rcpp -> R
    gsl_vector_view mult = gsl_matrix_row(myTest.smryStat, 0);
    gsl_vector_view Pmul = gsl_matrix_row(myTest.Psmry, 0);
    RcppGSL::VectorView umult = gsl_vector_subvector(&mult.vector, 1, nVars);
    RcppGSL::VectorView uPmul = gsl_vector_subvector(&Pmul.vector, 1, nVars);

    gsl_vector_view stat = gsl_matrix_column(myTest.smryStat, 0);
    gsl_vector_view Pstat = gsl_matrix_column(myTest.Psmry, 0);
    RcppGSL::VectorView signi = gsl_vector_subvector(&stat.vector, 1, nParam);
    RcppGSL::VectorView Psign = gsl_vector_subvector(&Pstat.vector, 1, nParam);

    RcppGSL::MatrixView usig = gsl_matrix_submatrix(myTest.smryStat, 1,1,nParam,nVars);
    RcppGSL::MatrixView uPsig= gsl_matrix_submatrix(myTest.Psmry, 1,1,nParam,nVars);

    List rs = List::create(
	Named("multstat" ) = gsl_vector_get(&mult.vector, 0),
	Named("Pmultstat") = gsl_vector_get(&Pmul.vector, 0),
        Named("unitmult" ) = umult,
        Named("Punitmult") = uPmul,
	Named("signific" ) = signi,
	Named("Psignific") = Psign,
	Named("unitsign" ) = usig,
	Named("Punitsign") = uPsig,
	Named("nSamp"    ) = myTest.nSamp,
	Named("aic"      ) = NumericVector(myTest.aic, myTest.aic+nVars)
    );

    // clear objects
    glmPtr[mtype]->releaseGlm();
    myTest.releaseTest();    
    gsl_vector_free(tm.smry_lambda);
    
    return rs;
}


// [[Rcpp::export]]
List RtoSmryCpp(const List & rparam,
                RcppGSL::Matrix & Y,
                RcppGSL::Matrix & X,
                Rcpp::Nullable<RcppGSL::Matrix> & bID)
{
    // pass regression parameters
    mv_Method mm;	
    mm.nboot = as<unsigned int>(rparam["nboot"]);
    mm.corr = as<unsigned int>(rparam["cor_type"]);
    mm.shrink_param = as<double>(rparam["shrink_param"]);
    mm.test = as<unsigned int>(rparam["test_type"]);
    mm.resamp = as<unsigned int>(rparam["resamp"]);
    mm.reprand = as<unsigned int>(rparam["reprand"]);
    mm.student = as<unsigned int>(rparam["studentize"]);
    mm.punit = as<unsigned int>(rparam["punit"]);
    mm.rsquare = as<unsigned int>(rparam["rsquare"]);

    unsigned int nRows = Y.nrow();
    unsigned int nVars = Y.ncol();
    unsigned int nParam = X.ncol();

    // do stuff	
//  clock_t clk_start, clk_end;
//  clk_start = clock();

    // initialize summary class
    Summary smry(&mm, Y, X);
	
    // Resampling indices
    if ( bID.isNotNull() ) {
        RcppGSL::Matrix m(bID);
        mm.nboot = m.nrow();	   
        smry.bootID = m;
    }    
  
// resampling test
    smry.resampTest();
//    smry.display();

    // Wrap gsl vectors with Rcpp, then Rcpp -> R 
    RcppGSL::VectorView umulti= gsl_matrix_row(smry.unitstat,0);
    RcppGSL::VectorView uPmult= gsl_matrix_row(smry.Punitstat,0);
    RcppGSL::MatrixView usigni= gsl_matrix_submatrix(smry.unitstat,1,0,nParam,nVars);
    RcppGSL::MatrixView uPsign= gsl_matrix_submatrix(smry.Punitstat,1,0,nParam,nVars);

    List rs = List::create(
	Named("multstat" )= smry.multstat[0],
	Named("Pmultstat")= smry.Pmultstat[0],
	Named("signific" )= NumericVector(smry.multstat+1, smry.multstat+nParam+1),
	Named("Psignific")= NumericVector(smry.Pmultstat+1, smry.Pmultstat+nParam+1),
        Named("unitmult" )= umulti,
        Named("Punitmult")= uPmult,
	Named("unitsign" )= usigni,
	Named("Punitsign")= uPsign,
	Named("nSamp"    )= smry.nSamp,
	Named("R2"       )= smry.R2
    );

    // clear objects
    smry.releaseSummary(); 
    
    return rs;
}

