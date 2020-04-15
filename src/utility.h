#ifndef UTILITY_H_
#define UTILITY_H_
#include <iostream>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h> 
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_matrix.h>

int my_spmatrix_scale_row_ccs_double(gsl_spmatrix * m, gsl_vector * x);

class BetaSpmSol {
  public:
	  BetaSpmSol(int nRows, int nParams, gsl_spmatrix* Xsp, gsl_spmatrix* XspT) {
		  WXsp_ = gsl_spmatrix_alloc(nRows,nParams);
		  Wsp_ = gsl_spmatrix_alloc(nRows, nRows);
		  XwXsp_ = gsl_spmatrix_alloc(nParams,nParams);
		  XrefspC_  = gsl_spmatrix_ccs(Xsp);
		  XrefspTC_ = gsl_spmatrix_ccs(XspT);
		  W_ = gsl_vector_alloc(nRows);
      const gsl_splinalg_itersolve_type *T_ = gsl_splinalg_itersolve_gmres;
      work_ = gsl_splinalg_itersolve_alloc(T_, nParams, 0);
      WXspC_ = gsl_spmatrix_ccs(WXsp_);
      XwXspC_ = gsl_spmatrix_ccs(XwXsp_);
		}
		void updateXspCCS(gsl_spmatrix* Xsp, gsl_spmatrix* XspT, int nParams) {
		  if (Xsp != NULL) {
				if (XrefspC_ != NULL) gsl_spmatrix_free(XrefspC_);
				XrefspC_  = gsl_spmatrix_ccs(Xsp);
			}
		  if (XspT != NULL) {
				if (XrefspC_ != NULL) gsl_spmatrix_free(XrefspTC_);
				XrefspTC_ = gsl_spmatrix_ccs(XspT);
			}
			gsl_splinalg_itersolve_free(work_);
      const gsl_splinalg_itersolve_type *T_ = gsl_splinalg_itersolve_gmres;
      work_ = gsl_splinalg_itersolve_alloc(T_, nParams, 0);
		}

		~BetaSpmSol() {
			gsl_spmatrix_free(WXsp_); 		
			gsl_spmatrix_free(Wsp_);
			gsl_spmatrix_free(XwXsp_);
			gsl_spmatrix_free(XrefspC_);
			gsl_spmatrix_free(XrefspTC_);
			gsl_vector_free(W_);
			gsl_splinalg_itersolve_free(work_);
		}

  public:
	  gsl_spmatrix* WXsp_; 		
		gsl_spmatrix* Wsp_;
    gsl_spmatrix* XwXsp_;
    gsl_spmatrix* XrefspC_;
    gsl_spmatrix* XrefspTC_;
	  gsl_vector* W_;
    gsl_splinalg_itersolve *work_;
    gsl_spmatrix* WXspC_;
    gsl_spmatrix* XwXspC_;
};

class MatrixSol {
  public:
	  MatrixSol(int nRows, int nParams) {
		  WX_ = gsl_matrix_alloc(nRows,nParams);
		  W_ = gsl_matrix_alloc(nRows, nRows);
		  XwX_ = gsl_matrix_alloc(nParams,nParams);
		  Xwz_ = gsl_vector_alloc(nParams);
		}
		~MatrixSol() {
			gsl_matrix_free(WX_);	
			gsl_matrix_free(W_);
			gsl_matrix_free(XwX_);
			gsl_vector_free(Xwz_);
		}
  public:
	  gsl_matrix* WX_; 		
		gsl_matrix* W_;
    gsl_matrix* XwX_;
	  gsl_vector* Xwz_;
};

class DisplayMatrix {
 public:	
	DisplayMatrix(gsl_matrix* m, const std::string& name) : mat_(m), name_(name) {}
	void Print() {
		int row = mat_->size1;
		int col = mat_->size2;
		std::cout << "matrix: " << name_ << "\n";
		for (int i = 0; i < row; ++i) {
			std::cout <<"[ ";
			for (int j = 0; j < col; ++j) {
			  std::cout << gsl_matrix_get(mat_, i, j) <<" ";
			}	
			std::cout << "]\n";
		}	
		std::cout<<"end of matrix "<< name_ <<"\n";
	}
	
	private:
		gsl_matrix* mat_;
		std::string name_;
};



#endif //end of unility_h_
