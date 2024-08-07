#include "abc_000_macro.h"
#include "abc_000_warning.h"

#include <math.h>
#include <string.h>

#include "abc_mcmc.h"
#include "abc_ts_func.h"
#include "abc_mat.h"
#include "abc_math.h"    //sumlog, sum_log_diagv2, fastlog
#include "abc_rand.h"
#include "abc_mem.h"
#include "abc_blas_lapack_lib.h"

#include "globalvars.h"
#include "beastv2_header.h"
#include "beastv2_prior_precfunc.h" 

#define MODEL (*model)

static void SetNtermsPerPrecGrp_prec01(I16PTR nTermsPerPrecGrp, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState) {
	// [0:ConstPrec, 1:UniformPrec, 2:ComponentWise], 3:OrderWise
	return;
}
static void SetNtermsPerPrecGrp_prec2(I16PTR nTermsPerPrecGrp, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState) {
	// [0:ConstPrec, 1:UniformPrec, 2:ComponentWise], 3:OrderWise
	for (int i = 0; i <  NUMBASIS; i++) {
		nTermsPerPrecGrp[i] = b[i].K;
	}
	return;
}
static void SetNtermsPerPrecGrp_prec3(I16PTR nTermsPerPrecGrp, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState) {
	// 0:ConstPrec, 1:UniformPrec, 2:ComponentWise, [3:OrderWise]
	memset(nTermsPerPrecGrp, 0,  precState->nPrecGrp *sizeof(I16));

	for (int id = 0; id < NUMBASIS; ++id) {
		I16PTR   _nTermsPerGrp = nTermsPerPrecGrp + b[id].offsetPrec;
		U08PTR   termType      = b[id].termType;
		for (int k = 0; k < b[id].K; ++k) {	 
			_nTermsPerGrp[ termType[k]-1 ]++;
		}
	}
}

// [0:ConstPrec, 1:UniformPrec], 2:ComponentWise, 3:OrderWise
void SetPrecXtXDiag_prec01(F32PTR precXtXDiag, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState) {return;}


void SetPrecXtXDiag_prec2(F32PTR precXtXDiag, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState) {
	// 0:ConstPrec, 1:UniformPrec, [2:ComponentWise], 3:OrderWise

	for (int id = 0; id < NUMBASIS; ++id) {
		//F32PTR    precVec = MODEL.precState.precVec + MODEL.b[id].offsetPrec;
		// For componentwise precision type, precVec takes only one value for each componet,
		// so id can be used as an index rather than using b[id].offsetPrec
		F32	 prec  = precState->precVec[id];  
		for (int k = 0; k < b[id].K; ++k) {
			*precXtXDiag++ = prec;
		}
	}
}

void SetPrecXtXDiag_prec3(F32PTR precXtXDiag, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState) {
	// 0:ConstPrec, 1:UniformPrec, 2:ComponentWise, [3:OrderWise]
	for (int id = 0; id <  NUMBASIS; ++id) {
		U08PTR    termType =  b[id].termType;
		F32PTR    prec     = precState->precVec + b[id].offsetPrec;
		for (int k = 0; k < b[id].K; ++k) {
			*precXtXDiag++ = prec[termType[k] - 1];
		}
	} 
}


void SetPropPrecXtXDiag_NtermsPerGrp_prec01(BEAST2_MODEL_PTR model, BEAST2_BASIS_PTR basis, NEWTERM_PTR new  ){
	// [0:ConstPrec, 1:UniformPrec], 2:ComponentWise, 3:OrderWise
	return; 
}

void SetPropPrecXtXDiag_NtermsPerGrp_prec2(BEAST2_MODEL_PTR model, BEAST2_BASIS_PTR basis, NEWTERM_PTR new ) {
	// 0:ConstPrec, 1:UniformPrec, [2:ComponentWise], 3:OrderWise

	// new is not used

	const int k1     = new->newcols.k1;
	const int k2_old = new->newcols.k2_old;
	const int k2_new = new->newcols.k2_new;

	const I32 Kold = MODEL.curr.K;
	SCPY(  k1 - 1,      MODEL.curr.precXtXDiag,           MODEL.prop.precXtXDiag);
	SCPY(Kold- k2_old,  MODEL.curr.precXtXDiag + k2_old,  MODEL.prop.precXtXDiag + k2_new);
	
	//F32PTR    precVec  = MODEL.precState.precVec     + basis->offsetPrec;
	//F32		prec	 = *precVec;

	// basis id is used to index the prec value
	const I32 id   = basis - model->b;
	const F32 prec = MODEL.precState.precVec[id];

	F32PTR  precXtXDiag = MODEL.prop.precXtXDiag;
	for (int k =  k1 ; k <= k2_new; ++k) {
		precXtXDiag[k - 1] = prec;		
	}

	// CHange the number of terms in the selected component
	// NUMBASIS=model.precState.nPrecGrp
	memcpy(MODEL.prop.nTermsPerPrecGrp, MODEL.curr.nTermsPerPrecGrp, sizeof(I16) * model->NUMBASIS);
	MODEL.prop.nTermsPerPrecGrp[id] +=  k2_new - k2_old;
}

void SetPropPrecXtXDiag_NtermsPerGrp_prec3(BEAST2_MODEL_PTR model, BEAST2_BASIS_PTR basis, NEWTERM_PTR new )
{


	int k1     = new->newcols.k1;
	int k2_old = new->newcols.k2_old;
	int k2_new = new->newcols.k2_new;

	 //  3:OrderWise

	I32 Kold = MODEL.curr.K;
	SCPY(k1   - 1,      MODEL.curr.precXtXDiag,          MODEL.prop.precXtXDiag);
	SCPY(Kold - k2_old, MODEL.curr.precXtXDiag + k2_old, MODEL.prop.precXtXDiag + k2_new);

	F32PTR  precXtXDiag = MODEL.prop.precXtXDiag  + k1 - 1;
	F32PTR  prec        = MODEL.precState.precVec + basis->offsetPrec;

	if (basis->type == SEASONID){
		for (int i = 0; i < new->numSeg; i++) {
			for (int order = new->SEG[i].ORDER1; order <= new->SEG[i].ORDER2; order++) {
				*precXtXDiag++ = prec[order - 1];
				*precXtXDiag++ = prec[order - 1];
			}
		}
	}
	else if (basis->type == TRENDID) {
		for (int i = 0; i < new->numSeg; i++) {
			for (int order = new->SEG[i].ORDER1; order <= new->SEG[i].ORDER2; order++) {
				*precXtXDiag++ = prec[order];
			}
		}
	}
	else if (basis->type == OUTLIERID) {
		//new->numSeg == Knewterm
		for (int i = 0; i < new->numSeg; i++)
			precXtXDiag[i] = prec[0];
	}

	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////

	#undef  NEW   //undo #define NEW(class_def)  R_do_new_object(class_def)
	#define NEW (*new)

	I16PTR nTermsPerPrecGrp = MODEL.prop.nTermsPerPrecGrp + basis->offsetPrec;
	memcpy(MODEL.prop.nTermsPerPrecGrp, MODEL.curr.nTermsPerPrecGrp, sizeof(I16) * MODEL.precState.nPrecGrp);
	memset(nTermsPerPrecGrp, 0, sizeof(I16) * basis->nPrec);


	if (basis->type == SEASONID || basis->type == TRENDID) {
		I32 k1          = k1     - basis->Kbase;
		I32 k2old       = k2_old - basis->Kbase;
		I32 k2new       = k2_new - basis->Kbase;
		U08PTR termType = basis->termType;

		for (int i = 1; i <= k1 - 1; i++) {
			nTermsPerPrecGrp[termType[i - 1] - 1]++;
		}

		if (basis->type == SEASONID) {			
			for (int i = 0; i < NEW.numSeg; i++) {
				for (int order = new->SEG[i].ORDER1; order <= new->SEG[i].ORDER2; order++) {
					nTermsPerPrecGrp[order - 1] += 2;
				}
			}
		} else		{
			for (int i = 0; i < NEW.numSeg; i++) {
				for (int order = new->SEG[i].ORDER1; order <= new->SEG[i].ORDER2; order++) {
					nTermsPerPrecGrp[order ] += 1;
				}
			}

		}

		I32 Kbasis_old = basis->K;
		for (int i = k2old + 1; i <= Kbasis_old; i++) {
			nTermsPerPrecGrp[termType[i - 1] - 1]++;
		}
	 
	}
	else if (basis->type == OUTLIERID) {
		nTermsPerPrecGrp[0] = NEW.nKnot_new;
	}
#undef  NEW

}

void ResamplePrecValues_prec0(BEAST2_MODEL_PTR model, BEAST2_HyperPar* hyperPar, VOID_PTR stream) {
	// 0:ConstPrec
	return;
}

void ResamplePrecValues_prec1(BEAST2_MODEL_PTR model, BEAST2_HyperPar * hyperPar, VOID_PTR stream) {
	// 1:UniformPrec

	I32 K    = MODEL.curr.K;
	F32 sumq = DOT(K, MODEL.beta, MODEL.beta);

	r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, (*(VSLStreamStatePtr*)stream), 1L, MODEL.precState.precVec, (hyperPar->del_1 + K * 0.5f), 0, 1.f);
	F32 newPrecVal                = MODEL.precState.precVec[0] / (hyperPar->del_2 + 0.5f * sumq / model->sig2[0]);
	MODEL.precState.precVec[0]    = newPrecVal > MIN_PREC_VALUE? newPrecVal: MIN_PREC_VALUE;
	MODEL.precState.logPrecVec[0] = logf(MODEL.precState.precVec[0]);
		 
}
void ResamplePrecValues_prec2(BEAST2_MODEL_PTR model, BEAST2_HyperPar * hyperPar, VOID_PTR stream) {
	
	// 0:ConstPrec, 1:UniformPrec, [2:ComponentWise], 3:OrderWise

		for (int id = 0; id < MODEL.NUMBASIS; id++) {
			BEAST2_BASIS_PTR	basis	= &MODEL.b[id];
			F32PTR				beta	= MODEL.beta + basis->Kbase;
			I32                 K		= basis->K;

			if (K <= 0) continue;

			F32		sumq = DOT(K, beta, beta);
			r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, (*(VSLStreamStatePtr*)stream), 1, &(MODEL.precState.precVec[id]), (hyperPar->del_1 + K * 0.5f), 0, 1.f);
			F32 newPrecVal          = MODEL.precState.precVec[id] / (hyperPar->del_2 + 0.5f * sumq / MODEL.sig2[0]);
			MODEL.precState.precVec[id]		= newPrecVal > MIN_PREC_VALUE ? newPrecVal : MIN_PREC_VALUE;
			MODEL.precState.logPrecVec[id]	= logf(MODEL.precState.precVec[id]);	
		}		
}

void ResamplePrecValues_prec3(BEAST2_MODEL_PTR model, BEAST2_HyperPar * hyperPar, VOID_PTR stream) {

	// 0:ConstPrec, 1:UniformPrec, 2:ComponentWise, [3:OrderWise]

	for (int id = 0; id < MODEL.NUMBASIS; id++) {
		BEAST2_BASIS_PTR basis = &MODEL.b[id];
		F32PTR prec     = MODEL.precState.precVec    + basis->offsetPrec;;
		F32PTR logPrec  = MODEL.precState.logPrecVec + basis->offsetPrec;;
		U08PTR termType = basis->termType;

		F32PTR beta	= MODEL.beta + basis->Kbase;
		for (int i = 1; i <= basis->nPrec; i++) {
			F32 sumq	= 0;
			I32 K		= 0;
			for (int j = 0; j < basis->K; j++) {
				if (termType[j] == i) {
					sumq += beta[j]*beta[j]; K++;
				}
			}
			if (K > 0) {
				r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, (*(VSLStreamStatePtr*)stream), 1, &prec[i - 1], (hyperPar->del_1 + K * 0.5f), 0, 1.f);
				F32 newPrecVal = prec[i - 1] / (hyperPar->del_2 + 0.5f * sumq / MODEL.sig2[0]);
				prec[i - 1]    = newPrecVal > MIN_PREC_VALUE ? newPrecVal : MIN_PREC_VALUE;
				logPrec[i - 1] = logf(prec[i - 1]);
			}
		}

	}
}

void IncreasePrecValues_prec0(BEAST2_MODEL_PTR model) {
	// 0:ConstPrec
	return;
}

void IncreasePrecValues_prec1(BEAST2_MODEL_PTR model) {
	// 1:UniformPrec
	MODEL.precState.precVec[0]    =  2 * MODEL.precState.precVec[0] ;
	MODEL.precState.logPrecVec[0] = logf(MODEL.precState.precVec[0]);

}
 
void IncreasePrecValues_prec2(BEAST2_MODEL_PTR model) {
	// 0:ConstPrec, 1:UniformPrec, [2:ComponentWise], 3:OrderWise
	for (int id = 0; id < MODEL.NUMBASIS; ++id) {
		MODEL.precState.precVec[id]    = MODEL.precState.precVec[id]+ MODEL.precState.precVec[id];
		MODEL.precState.logPrecVec[id] = logf(MODEL.precState.precVec[id]);
	}

}

void IncreasePrecValues_prec3(BEAST2_MODEL_PTR model) {

	// 0:ConstPrec, 1:UniformPrec, 2:ComponentWise, [3:OrderWise]

	for (int id = 0; id < MODEL.NUMBASIS; ++id) {
		BEAST2_BASIS_PTR basis = &MODEL.b[id];
		F32PTR prec     = MODEL.precState.precVec    + basis->offsetPrec;;
		F32PTR logPrec  = MODEL.precState.logPrecVec + basis->offsetPrec;;
		U08PTR termType = basis->termType;

		for (int i = 1; i <= basis->nPrec; i++) {
			I32 K   = 0;
			for (int j = 0; j < basis->K; j++) {
				if (termType[j] == i) { K++; }
			}
			if (K > 0) {								
				prec[i - 1]    = prec[i - 1] + prec[i - 1];
				logPrec[i - 1] = logf(prec[i - 1]);
			}
		}

	}
}


void ComputeMargLik_prec01(BEAST2_MODELDATA_PTR data, PRECSTATE_PTR precState,	BEAST2_YINFO_PTR yInfo, 
								 BEAST2_HyperPar_PTR hyperPar)
{	
	 I32 K = data->K;
	 solve_U_as_LU_invdiag_sqrmat(data->cholXtX, data->XtY, data->beta_mean, K);
	//Compute alpha2_star: Two ways to compute alpha2_star: YtY-m*V*m
	/* ---THE FIRST WAY---
			//GlobalMEMBuf_1st = Xnewterm + K_newTerm*Npad;
			r_cblas_scopy(KNEW, beta_mean_prop, 1, GlobalMEMBuf_2nd, 1);
			//cblas_dtrmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
			r_cblas_strmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, KNEW, cholXtX_prop, KNEW, GlobalMEMBuf_2nd, 1);
			//double cblas_ddot (const MKL_INT n, const double *x, const MKL_INT incx, const double *y, const MKL_INT incy);
			basis_prop->alpha2_star = yInfo.YtY - DOT(KNEW, GlobalMEMBuf_2nd, GlobalMEMBuf_2nd);
	*/
	/*---THE SECOND WAY-- */
	F32 alpha2_star		 =  (yInfo->YtY_plus_alpha2Q[0] - DOT(K, data->XtY,data->beta_mean))*0.5;

	//half_log_det_post  = sum(log(1. / diag(U)))	 
	//half_log_det_post  = -sum_log_diagv2(MODEL.prop.cholXtX, KNEW);
	F32 half_log_det_post = sum_log_diagv2(data->cholXtX, K);

	////Copy the diagonal of P_U to buf
	// vmsLn(KNEW, GlobalMEMBuf_1st, GlobalMEMBuf_1st, VML_HA);
	// ippsSum_32f(GlobalMEMBuf_1st, KNEW, &FLOAT_SHARE.half_log_det_post, ippAlgHintFast);
	// r_ippsSumLn_32f(GlobalMEMBuf_2nd, KNEW, &half_log_det_post);	

	//half_log_det_prior = -0.5 * (k_SN * log(prec(1)) + length(k_const_terms) * log(prec(2)) + length(linear_terms) * log(prec(3)));
	// F32 half_log_det_prior = -0.5f * KNEW * *modelPar.LOG_PREC[0]

	F32 half_log_det_prior = -0.5f* precState->logPrecVec[0] * K;

	//log_ML = half_log_det_post - half_log_det_prior - (n / 2 + b) * log(a + a_star / 2);
	F32 marg_lik   = half_log_det_post  -half_log_det_prior	- yInfo->alpha1_star * fastlog(alpha2_star);

	data->alpha2Q_star[0] = alpha2_star;
	data->marg_lik        = marg_lik;

	//if (q == 1) {// added for MRBEAST
	data->alpha2Q_star[0] = max(data->alpha2Q_star[0], MIN_ALPHA2_VALUE);
	//}
}
 
void ComputeMargLik_prec23(BEAST2_MODELDATA_PTR data, PRECSTATE_PTR precState,BEAST2_YINFO_PTR yInfo,
	BEAST2_HyperPar_PTR hyperPar) {

	I32 K = data->K;
	solve_U_as_LU_invdiag_sqrmat(data->cholXtX, data->XtY, data->beta_mean, K);
	//Compute alpha2_star: Two ways to compute alpha2_star: YtY-m*V*m
	/* ---THE FIRST WAY---
			//GlobalMEMBuf_1st = Xnewterm + K_newTerm*Npad;
			r_cblas_scopy(KNEW, beta_mean_prop, 1, GlobalMEMBuf_2nd, 1);
			//cblas_dtrmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
			r_cblas_strmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, KNEW, cholXtX_prop, KNEW, GlobalMEMBuf_2nd, 1);
			//double cblas_ddot (const MKL_INT n, const double *x, const MKL_INT incx, const double *y, const MKL_INT incy);
			basis_prop->alpha2_star = yInfo.YtY - DOT(KNEW, GlobalMEMBuf_2nd, GlobalMEMBuf_2nd);
	*/

	/*---THE SECOND WAY-- */
	F32 alpha2_star        = (yInfo->YtY_plus_alpha2Q[0] - DOT(K, data->XtY, data->beta_mean)) * 0.5;

	//half_log_det_post  = sum(log(1. / diag(U)))	 
	//half_log_det_post  = -sum_log_diagv2(MODEL.prop.cholXtX, KNEW);
	F32 half_log_det_post = sum_log_diagv2(data->cholXtX, K);

	////Copy the diagonal of P_U to buf
	// vmsLn(KNEW, GlobalMEMBuf_1st, GlobalMEMBuf_1st, VML_HA);
	// ippsSum_32f(GlobalMEMBuf_1st, KNEW, &FLOAT_SHARE.half_log_det_post, ippAlgHintFast);
	// r_ippsSumLn_32f(GlobalMEMBuf_2nd, KNEW, &half_log_det_post);	

	//half_log_det_prior = -0.5 * (k_SN * log(prec(1)) + length(k_const_terms) * log(prec(2)) + length(linear_terms) * log(prec(3)));
	// F32 half_log_det_prior = -0.5f * KNEW * *modelPar.LOG_PREC[0]

	/***********************PrecPiorType=2*******************************/	 
	F32  half_log_det_prior = 0;
	for (I32 i = 0; i < precState->nPrecGrp; i++)
		half_log_det_prior += precState->logPrecVec[i] * data->nTermsPerPrecGrp[i];
	half_log_det_prior *= -0.5;
	/***********************PrecPiorType=2*******************************/

	//log_ML = half_log_det_post - half_log_det_prior - (n / 2 + b) * log(a + a_star / 2);
	F32 marg_lik = half_log_det_post - half_log_det_prior -yInfo->alpha1_star * fastlog(alpha2_star);

	data->alpha2Q_star[0] = alpha2_star;
	data->marg_lik        = marg_lik;

	//if (q == 1) {// added for MRBEAST
	data->alpha2Q_star[0] = max(data->alpha2Q_star[0], MIN_ALPHA2_VALUE);
	//}
}

 
void MR_ComputeMargLik_prec01(BEAST2_MODELDATA_PTR data, PRECSTATE_PTR precState,  BEAST2_YINFO_PTR yInfo,
								BEAST2_HyperPar_PTR hyperPard)
{
	I32 q = yInfo->q;
	I32 K = data->K;

	//MRBEAST
	F32PTR beta_mean = data->beta_mean;
	F32PTR XtY       = data->XtY;
	//BEAST: solve_U_as_LU_invdiag_sqrmat(cholXtX, XtY, beta_mean, K);
	//MRBEAST--------start
	solve_U_as_LU_invdiag_sqrmat_multicols(data->cholXtX, XtY, beta_mean, K, q);
	//MRBEAST--------end
	
	//BEAST: F32 alpha2_star = (pyInfo->YtY_plus_alpha2Q[0] - DOT(K,  XtY,  beta_mean)) * 0.5;
    //MRBEAST--------start
	//cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, k_newTerm, 1, N, 1, X_mars_prop + (k1_new - 1)*N, N, Y, N, 0, buff1, k_newTerm);
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, K, 1.f, beta_mean, K, XtY, K, 0.f, data->alpha2Q_star, q);
	r_ippsSub_32f(data->alpha2Q_star, yInfo->YtY_plus_alpha2Q, data->alpha2Q_star, q * q);	
	//MRBEAST--------end
	

	//half_log_det_post; = sum(log(1. / diag(U)))
	F32 half_log_det_post = sum_log_diagv2(data->cholXtX, K);

	 //half_log_det_prior = -0.5 * (k_SN * log(prec(1)) + length(k_const_terms) * log(prec(2)) + length(linear_terms) * log(prec(3)));	
	F32 half_log_det_prior = -0.5f * K * precState->logPrecVec[0];
 
	/****************************************************************************************************************************/
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', q, data->alpha2Q_star, q); // Choleskey decomposition; only the upper triagnle elements are used

	//F32 log_det_alphaQ = sum_log_diagv2(curmodel->alpha2Q_star, q);
	F32 sumLn_alphaQ_det = sum_log_diagv2(data->alpha2Q_star, q);
 
	//log_ML = half_log_det_post - half_log_det_prior - (n / 2 + b) * log(a + a_star / 2);
	data->marg_lik = q*(half_log_det_post - half_log_det_prior) - yInfo->alpha1_star * sumLn_alphaQ_det*2.f; 	 

	//r_printf("Lik: det_post %f\n", f32_abs_sum(beta_mean, K * q));
	//r_printf("Lik: det_post %f\n", f32_sum_matrixdiag(data->cholXtX, K));
}



void SetUpPrecFunctions(I08 precPriorType, I32 q, PREC_FUNCS * funcs) {

	
	if (q == 1) {// for BEASTv4 
			if (  precPriorType == ConstPrec || precPriorType == UniformPrec) {
				funcs->SetPrecXtXDiag              = SetPrecXtXDiag_prec01;
				funcs->SetNtermsPerPrecGrp         = SetNtermsPerPrecGrp_prec01;
				funcs->SetPropPrecXtXDiag_NtermsPerGrp = SetPropPrecXtXDiag_NtermsPerGrp_prec01;
				funcs->chol_addCol                     = chol_addCol_skipleadingzeros_prec_invdiag;
			}
			else if (precPriorType == ComponentWise)
			{
				funcs->SetPrecXtXDiag                  = SetPrecXtXDiag_prec2;
				funcs->SetNtermsPerPrecGrp             = SetNtermsPerPrecGrp_prec2;
				funcs->SetPropPrecXtXDiag_NtermsPerGrp = SetPropPrecXtXDiag_NtermsPerGrp_prec2;
				funcs->chol_addCol                     = chol_addCol_skipleadingzeros_precVec_invdiag;

			}
			else if (precPriorType == OrderWise )
			{

				funcs->SetPrecXtXDiag                  = SetPrecXtXDiag_prec3;
				funcs->SetNtermsPerPrecGrp             = SetNtermsPerPrecGrp_prec3;
				funcs->SetPropPrecXtXDiag_NtermsPerGrp = SetPropPrecXtXDiag_NtermsPerGrp_prec3;
				funcs->chol_addCol                     = chol_addCol_skipleadingzeros_precVec_invdiag;
		 
			}

			if      (precPriorType == 0)	funcs->IncreasePrecValues = IncreasePrecValues_prec0;
			else if (precPriorType == 1) 	funcs->IncreasePrecValues = IncreasePrecValues_prec1;
			else if (precPriorType == 2)	funcs->IncreasePrecValues = IncreasePrecValues_prec2;
			else if (precPriorType == 3)	funcs->IncreasePrecValues = IncreasePrecValues_prec3;

			if (precPriorType == 0)			funcs->ResamplePrecValues = ResamplePrecValues_prec0;
			else if (precPriorType == 1) 	funcs->ResamplePrecValues = ResamplePrecValues_prec1;
			else if (precPriorType == 2)	funcs->ResamplePrecValues = ResamplePrecValues_prec2;
			else if (precPriorType == 3)	funcs->ResamplePrecValues = ResamplePrecValues_prec3;

			if (precPriorType == 0)	        funcs->ComputeMargLik = ComputeMargLik_prec01;
			else if (precPriorType == 1) 	funcs->ComputeMargLik = ComputeMargLik_prec01;
			else if (precPriorType == 2)	funcs->ComputeMargLik = ComputeMargLik_prec23;
			else if (precPriorType == 3)	funcs->ComputeMargLik = ComputeMargLik_prec23;


	}

	if (q> 1) {// for BEASTv4 
			if (  precPriorType == ConstPrec || precPriorType == UniformPrec) {
				funcs->SetPrecXtXDiag              = SetPrecXtXDiag_prec01;
				funcs-> SetNtermsPerPrecGrp       = SetNtermsPerPrecGrp_prec01;
				funcs->SetPropPrecXtXDiag_NtermsPerGrp = SetPropPrecXtXDiag_NtermsPerGrp_prec01;
				funcs->chol_addCol                     = chol_addCol_skipleadingzeros_prec_invdiag;
			}
			else if (precPriorType == ComponentWise)
			{
				funcs->SetPrecXtXDiag              = SetPrecXtXDiag_prec2;
				funcs-> SetNtermsPerPrecGrp       = SetNtermsPerPrecGrp_prec2;
				funcs->SetPropPrecXtXDiag_NtermsPerGrp = SetPropPrecXtXDiag_NtermsPerGrp_prec2;
				funcs->chol_addCol                     = chol_addCol_skipleadingzeros_precVec_invdiag;

			}
			else if (precPriorType == OrderWise )
			{

				funcs->SetPrecXtXDiag              = SetPrecXtXDiag_prec3;
				funcs-> SetNtermsPerPrecGrp       = SetNtermsPerPrecGrp_prec3;
				funcs->SetPropPrecXtXDiag_NtermsPerGrp = SetPropPrecXtXDiag_NtermsPerGrp_prec3;
				funcs->chol_addCol                     = chol_addCol_skipleadingzeros_precVec_invdiag;
		 
			}

			if      (precPriorType == 0)	funcs->IncreasePrecValues = IncreasePrecValues_prec0;
			else if (precPriorType == 1) 	funcs->IncreasePrecValues = IncreasePrecValues_prec1;
			else if (precPriorType == 2)	funcs->IncreasePrecValues = IncreasePrecValues_prec2;
			else if (precPriorType == 3)	funcs->IncreasePrecValues = IncreasePrecValues_prec3;

			if (precPriorType == 0)			funcs->ResamplePrecValues = ResamplePrecValues_prec0;
			else if (precPriorType == 1) 	funcs->ResamplePrecValues = ResamplePrecValues_prec1;
			else if (precPriorType == 2)	funcs->ResamplePrecValues = ResamplePrecValues_prec2;
			else if (precPriorType == 3)	funcs->ResamplePrecValues = ResamplePrecValues_prec3;

			if (precPriorType == 0)	        funcs->ComputeMargLik = MR_ComputeMargLik_prec01;  // For MRBEAST
			else if (precPriorType == 1) 	funcs->ComputeMargLik = MR_ComputeMargLik_prec01;// ComputeMargLik_prec01;
			else if (precPriorType == 2)	funcs->ComputeMargLik = MR_ComputeMargLik_prec01; //ComputeMargLik_prec2;
			else if (precPriorType == 3)	funcs->ComputeMargLik = MR_ComputeMargLik_prec01; //ComputeMargLik_prec3;


	}
}

#include "abc_000_warning.h"
