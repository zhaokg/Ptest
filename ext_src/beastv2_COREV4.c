#include "abc_000_macro.h"
#include "abc_000_warning.h"

#if defined(COMPILER_MSVC)
#include "intrin.h"                //_rdstc
#endif

#include <string.h>	               //memset memcpy
#include <time.h>
#include <math.h>

#include "abc_001_config.h"
#include "abc_mem.h"              // Independent of R/Matlab,  VC/GNU, or MY/MKL.
#include "abc_blas_lapack_lib.h"  // Slight dependence on the choice of VC/GNU. Dependence on MY/MKL. Independacne of R/Matlab.
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_timer.h"
#include "abc_mcmc.h"
#include "abc_mat.h"
#include "abc_rand.h"
#include "abc_vec.h"   // for f32_add_v_v2_vec_in_place, f32_diff_back,i32_increment_bycon_inplace i32_to_f32_scaelby_inplace, f32_sx_sxx_toavstd_inplace 
#include "abc_math.h"  // for fastexp, fastsqrt only

#include <stdio.h>	               //fprintf fopen FILE #include<stdio.h>  // Need _GNU_SOURCE for manylinux; otherwise report /usr/include/stdio.h:316:6: error: unknown type name '_IO_cookie_io_functions_t'

#include "globalvars.h"  

#include "beastv2_header.h"
#include "beastv2_func.h" 
#include "beastv2_model_allocinit.h" 
#include "beastv2_prior_precfunc.h" 
#include "beastv2_xxyy_allocmem.h" 
#include "beastv2_io.h" 

//#include <unistd.h> // char* getcwd(char* buf, size_t size);

#define LOCAL(...) do{ __VA_ARGS__ } while(0);
//extern MemPointers* mem;
//time_t start, end; 

#define  DEBUG_MODE  0

int beast2_main_corev4(void)   {

	
	// A struct to track allocated pointers   
	// Do not use 'const MemPointers MEM' bcz Clang will asssume other fields as zeros (e.g., alloc, and alloc0).
	MemPointers MEM = (MemPointers){.init = mem_init,};
	MEM.init(&MEM);

#if DEBUG_MODE == 1
	MEM.checkHeader = 1;
	//mem = &MEM;
#endif

	VSLStreamStatePtr stream;

	// Get the Option parameters from the global pointer GLOBAL_OPTIONS, which was already 
	// set prior to this point (e.g., time-series length N, and number of pixels).	
	const BEAST2_OPTIONS_PTR	opt   = GLOBAL_OPTIONS;
	const BEAST2_EXTRA          extra = opt->extra;
	
	typedef int QINT;
	//const   QINT  q = 1L;
	const QINT  q   = opt->io.q;
	
	// Allocate MEMORY FOR BASIS VARIABLE: Initialzie two pointers to BASIS
	BEAST2_MODEL  MODEL = {0,};
	AllocInitModelMEM(&MODEL, opt, &MEM);

	//Initializing the random number generaotr	
	LOCAL( 	
		U64 seed = (opt->mcmc.seed == 0) ? TimerGetTickCount() : (opt->mcmc.seed+0x4f352a3dc);
		r_vslNewStream(&stream, VSL_BRNG_MT19937, seed);   	
	)

	const U32PTR  RND32        = MyALLOC(MEM, MAX_RAND_NUM,     U32, 64);
	const U16PTR  RND16        = MyALLOC(MEM, MAX_RAND_NUM * 2, U16, 64);
	const U08PTR  RND08        = MyALLOC(MEM, MAX_RAND_NUM * 4, U08, 64);	
	const F32PTR  RNDGAMMA     = MyALLOC(MEM, MAX_RAND_NUM,     F32, 64);	
	const U32PTR  RND32_END    = RND32		+ MAX_RAND_NUM - 7;
	const U16PTR  RND16_END    = RND16		+ MAX_RAND_NUM * 2 - 7;
	const U08PTR  RND08_END    = RND08		+ MAX_RAND_NUM * 4 - 7 -3;     //-3 bcz GenRandomBasis will also consume the stream bits
	const F32PTR  RNDGAMMA_END = RNDGAMMA	+ MAX_RAND_NUM - MODEL.precState.nPrecGrp-1L;
		
	// Allocate mem for current covariate/design matrix (Xt_mars), proposed new terms (Xnewterm),
	// and subset matrix corresponding to rows of missing values.		
	const F32PTR Xt_mars;
	const F32PTR Xnewterm;      //will be re-used as a temp mem for multiple purposes
	const F32PTR Xt_zeroBackup; //mem for saving Xrows of the missing rows	
	AllocateXXXMEM(&Xt_mars, &Xnewterm, &Xt_zeroBackup,&MODEL,opt,&MEM);

	// yInfo used to save the current time series to be processed
	BEAST2_YINFO     yInfo;
	AllocateYinfoMEM(&yInfo, opt, &MEM);
	
	// Allocate the output memory for a single chain (resultChain) and the averaged
	// result of all chains ( result)
	
	// Cannot make it constant; otehrwise resultChina will be treated as NULL and cause errors in mempcy(resultChina.xx)
	//beastv2_COREV4.c:677:13: warning: null passed to a callee that requires a non-null argument [-Wnonnull]
    BEAST2_RESULT resultChain = { NULL, };
	BEAST2_RESULT result      = { NULL, };
	BEAST2_Result_AllocMEM(&resultChain, opt, &MEM); 	
	BEAST2_Result_AllocMEM(&result,      opt, &MEM);
	


	// Pre-allocate memory to save samples for calculating credibile intervals	 
	const   I32  NumCIVars = MODEL.NUMBASIS + opt->extra.computeTrendSlope;
	 CI_PARAM     ciParam   = { 0, };
	CI_RESULT    ci[MAX_NUM_BASIS + 1];
	if (extra.computeCredible) {
		ConstructCIStruct(opt->mcmc.credIntervalAlphaLevel, opt->mcmc.samples, opt->io.N * opt->io.q,  //for MRBEAST
			              NumCIVars, &MEM, &extra.fastCIComputation, &ciParam, ci);
	}

	if (extra.computeCredible) {
		I32  Npad           = (opt->io.N + 7) / 8 * 8;
		I32  XnewtermOffset = 0;
		Npad = opt->io.N;    //Correct for the inconsitency of X and Y in gemm and gemv

		I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
		I08 hasTrendCmpnt   = 1;
		I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
		
		int numCIVars      =  0;
		if (hasSeasonCmpnt) {
			ci[numCIVars].result     = resultChain.sCI;		         //season		
			ci[numCIVars].newDataRow = Xnewterm + XnewtermOffset;	 //season		
			numCIVars++;
			XnewtermOffset += Npad * q;   //FOR MRBEAST
		}

		if (hasTrendCmpnt) {
			ci[numCIVars].result      = resultChain.tCI;               //trend			
			ci[numCIVars].newDataRow  = Xnewterm  + XnewtermOffset;    //trend			
			numCIVars++;
			XnewtermOffset += Npad * q;   //FOR MRBEAST
		}	 

		if (hasOutlierCmpnt) {
		  ci[numCIVars].result     = resultChain.oCI,               //outlier
		  ci[numCIVars].newDataRow = Xnewterm + XnewtermOffset;     //outlier  
		  numCIVars++;
		  XnewtermOffset += Npad * q;   //FOR MRBEAST
		}

		if (opt->extra.computeTrendSlope) {
			ci[numCIVars].result     = resultChain.tslpCI;           //trend  slope		
			ci[numCIVars].newDataRow = Xnewterm + XnewtermOffset;    //trend  slope	
			numCIVars++;
			XnewtermOffset += Npad * q;   //FOR MRBEAST
		}
		//NumCAIvars should equal  MODEL.NUMBASIS + opt->extra.computeTrendSlope;

	} //NUMVAR_FOR_CI=3

	const CORESULT coreResults[MAX_NUM_BASIS];
	SetupPointersForCoreResults(coreResults, MODEL.b, MODEL.NUMBASIS, &resultChain);
		 
	const BEAST2_HyperPar  hyperPar = { .alpha_1=opt->prior.alpha1,.alpha_2=opt->prior.alpha2,.del_1=opt->prior.delta1, .del_2=opt->prior.delta2};

	/****************************************************************************/
	//		THE OUTERMOST LOOP: Loop through all the time series one by one
	/****************************************************************************/
	// Get conversion factors from counts to seceonds
	InitTimerFunc();
	StartTimer();
	SetBreakPointForStartedTimer();

	const PREC_FUNCS precFunc;
	SetUpPrecFunctions(opt->prior.precPriorType, opt->io.q, &precFunc);

	// Print a blank line to be backspaced by the follow
	if (extra.printProgress) {
		F32 frac = 0.0; I32 firstTimeRun = 1;
		printProgress1(frac, extra.consoleWidth, Xnewterm, firstTimeRun);
	}

	//#define __DEBUG__
	//#undef  __DEBUG__ 

	#if DEBUG_MODE ==1
		// Allocate a mem block and memset it to zero
		I32    N          = opt->io.N;
		I32    Npad       = (N + 7) / 8 * 8; Npad =  N;//Correct for the inconsitency of X and Y in gemm and gemv
		F32PTR flagSat    = MyALLOC0(MEM, N, I32, 64);
		F32PTR Xdebug     = MyALLOC0(MEM, Npad*(opt->prior.K_MAX+ opt->prior.K_MAX), I32, 64); // one Kmax for Xt_mars, and another for Xtbackup
	#endif

	//#define XtX_ByGroup XtX_ByGroup_FULL
	//#define MatxMat     MatxMat_FULL
	//#define MatxVec     MatxVec_FULL

	// Calculate the total number of time sereies to be processed
	const U32  NUM_PIXELS    = opt->io.numOfPixels;
	const U32  MCMC_SAMPLES  = opt->mcmc.samples;
	const U32  MCMC_THINNING = opt->mcmc.thinningFactor;
	const U32  MCMC_BURNIN   = opt->mcmc.burnin;
	const U32  MCMC_CHAINNUM = opt->mcmc.chainNumber;
	const U16  SEASON_BTYPE  = opt->prior.seasonBasisFuncType;

	void (*Update_XtX_from_Xnewterm)(F32PTR X, F32PTR Xnewterm, F32PTR XtX, F32PTR XtXnew, NEWTERM * NEW, BEAST2_MODEL * MODEL);
	void (*Update_XtY_from_Xnewterm)(F32PTR Y, F32PTR Xnewterm, F32PTR XtY, F32PTR XtYnew, NEWTERM * new, I32 q);

	U16  GROUP_MatxMat = (MODEL.sid < 0  || opt->prior.seasonBasisFuncType  != 3 ) && (MODEL.vid < 0  || opt->prior.trendBasisFuncType   != 3 )
		              && (MODEL.tid < 0  || opt->prior.trendBasisFuncType   != 2 ) && (MODEL.oid < 0  || opt->prior.outlierBasisFuncType != 2 );
 
	Update_XtX_from_Xnewterm = GROUP_MatxMat ? Update_XtX_from_Xnewterm_ByGroup : Update_XtX_from_Xnewterm_NoGroup;
	Update_XtY_from_Xnewterm = GROUP_MatxMat ? Update_XtY_from_Xnewterm_ByGroup : Update_XtY_from_Xnewterm_NoGroup;

	NUM_OF_PROCESSED_GOOD_PIXELS  = 0;  // this is a global variable.
	NUM_OF_PROCESSED_PIXELS       = 0;  // this is also a global variable.

	for (U32 pixelIndex = 1; pixelIndex <= NUM_PIXELS; pixelIndex++)
	{
		// Fecth a new time-series: set up Y, nMissing,  n, rowsMissing		
		F32PTR MEMBUF           = Xnewterm; // Xnewterm is a temp mem buf.
		BEAST2_fetch_timeSeries(&yInfo, pixelIndex,  MEMBUF, &(opt->io));

		F32PTR  Xtmp             = Xt_mars;
		U08     skipCurrentPixel = BEAST2_preprocess_timeSeries(&yInfo, MODEL.b, Xtmp, opt);		

        #if DEBUG_MODE ==1
			I32 accS[5]  = { 0,},   accT[5]  = { 0, };
			I32 flagS[5] = { 0, },  flagT[5] = { 0, };
			for (int i = 0; i < yInfo.nMissing; i++) { flagSat[yInfo.rowsMissing[i]] = getNaN();}
        #endif

		#define __START_IF_NOT_SKIP_TIMESESIRIES__    
		#define __END_IF_NOT_SKIP_TIMESESIRIES__                        

		__START_IF_NOT_SKIP_TIMESESIRIES__  
		if (!skipCurrentPixel) {

		if (q == 1) {  // for BEASTV4

				// alpha2_star  = alpha_2 + 0.5(YtY-beta*X'Y) = 0.5* (  [YtY+2*alpha_2] - beta*X'Y  )
				// YtY+2*alpha_2  is pre-cacluated here. THe "2" before alpha_2 is to undo the division
				// later in the calcution of alpha2_star
				yInfo.YtY_plus_alpha2Q[0] = yInfo.YtY_plus_alpha2Q[0] + 2 * hyperPar.alpha_2;
				
				//Pre-compute alpha1_start, which depends only on n.
				yInfo.alpha1_star = yInfo.n * 0.5 + hyperPar.alpha_1;    

		}	else {	// For MRBEAST

				 //YtY was computed from MV_Fetch; here alpaha_Q is added to its diagonal.			
				f32_add_val_matrixdiag(yInfo.YtY_plus_alpha2Q, hyperPar.alpha_2, q);

				yInfo.alpha1_star = yInfo.n * 0.5 + (hyperPar.alpha_1 + q - 1) * 0.5;
		}

		/****************************************************************************/
		// GENERATE A STREAM OF RANDOM NUMBERS FOR FUTURE USE
		/****************************************************************************/
		BEAST2_RNDSTREAM  RND;
		{
			RND.rnd32 = RND32, RND.rnd16 = RND16, RND.rnd08 = RND08, RND.rndgamma = RNDGAMMA;
			//vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE,  stream, MAX_RAND_NUM, rnd, 0, 1.0);		
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND32);
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND16);
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND08);
			// rndgamma needed to resmaple Sig2, depending on n and alpha_1 
			// Not used to sample prec bcz the degree of freedom there depends on the number of terms Kterms
			r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, MAX_RAND_NUM, RNDGAMMA, ( hyperPar.alpha_1 + yInfo.n * 0.5f), 0, 1);
		}
	 
		if (extra.dumpMCMCSamples) {
			result.smcmc = opt->io.out.result->smcmc;
			result.tmcmc = opt->io.out.result->tmcmc;
			result.omcmc = opt->io.out.result->omcmc;
		}

		// Clear up and zero out RESULT for initialization	 
		BEAST2_Result_FillMEM(&result, opt, 0);		

		// Make sure the initial preValues in MODEL.precState.precVec contians no NANs as
		// it will contain residual values from the previous time series
		ReInit_PrecValues(&MODEL, opt);

		/****************************************************************************/
		//Iterate all the chains. The individual chain result will be saved into resultChain
		/****************************************************************************/
		for ( U32 chainNumber= 0;  chainNumber < MCMC_CHAINNUM; chainNumber++)
		{
			const I32  N      = opt->io.N; 
			const I32  Npad   = N;  (N + 7) / 8 * 8; //Correct for the inconsitency of X and Y in gemm and gemv
			const I32  Npad16 = (N + 15) /16 * 16;	
			/****************************************************************************/
			//                 GENERATE AN INITIAL MODEL
			/****************************************************************************/			
			{   
				// Generate random knots (numknot,KNOTs and ORDERS)				
				GenarateRandomBasis(MODEL.b, MODEL.NUMBASIS, N, &RND,&yInfo);//CHANGE: nunKnot, ORDER, KNOT, K, KBase, Ks, Ke, or TERM_TYPE				
				 
				/*
				for (int i = 0; i < MODEL.NUMBASIS; ++i) {
					MODEL.b[i].nKnot     =   0;
					MODEL.b[i].KNOT[0]   = N + 1;
					MODEL.b[i].ORDER[0]  = MODEL.b[i].prior.minOrder;
					MODEL.b[i].CalcBasisKsKeK_TermType(&MODEL.b[i]); 	/// Get Ks and Ke for individula segments of each components
				}
				*/
			
				//Evaluate the initial model and compute its marg lik. CHANGE: DERIVE XMARS, K, BETA, BETA_MEAN, MARG_LIK
				// Xtmars is a cotinguous mem block consiting of Xtmars, Xnewterm, and Xt_zerobackup. The first part will
				// be filled with Xtmars, and the rest will be used as a temp block in this function call.
				// We don't use Xnewterm  or Xt_zerobackup as a temp mem buf bcz the zeroOutXmars function may need a much
				// larger  mem due to the many terms of the inital random model	
				// basis->K won't be updated inside the function and the old values from CalcBasisKsKeK is kept
 
				BEAST2_EvaluateModel(&MODEL.curr, MODEL.b, Xt_mars, N, MODEL.NUMBASIS, &yInfo, &hyperPar, &MODEL.precState, &precFunc); // for both BEAST and MRBEAST
			 
		
			}
		

			{
				// Find candidate positions for SEASON AND TREND				
				// nMissing & rowsMissing used for the outlier function
				CvtKnotsToBinVec(MODEL.b, MODEL.NUMBASIS, N, &yInfo);

				// Reset all knots to ones bcz the real extrem positiosn won't be updated till samples > 1
				memset(MODEL.extremePosVec, 1, N);
				for (I32 i = 0; i < yInfo.nMissing; ++i) MODEL.extremePosVec[yInfo.rowsMissing[i]] = 0;
				MODEL.extremPosNum = yInfo.n;

				// Initialize the deviation vector to a large initial value: Used only for the OUTLIER proposal function
				// Deviation will be updated once sample > 1
				f32_fill_val(1e30, MODEL.deviation, N);
				for (I32 i = 0; i < yInfo.nMissing; ++i) MODEL.deviation[yInfo.rowsMissing[i]]   = getNaN();
				MODEL.avgDeviation[0] = 1.0;

				// Clear up and zero out resultChain for initialization
				BEAST2_Result_FillMEM(&resultChain, opt, 0);

				// Update the Kbase for the bases after the 1st one. The first one is fixed at Kbase=0.					
				// Kbases needed in GetInfoBandList (for Update_XtX), SetPropPrecXtXDiag_NtermsPerGrp_prec3, and computeY
				MODEL.b[0].Kbase = 0;                           // This is the first time and also the only time b[0].Kbase is specified
				UpdateBasisKbase(MODEL.b, MODEL.NUMBASIS, 0);	//CHANGE: MODEL.b[i].Kbase for i>basisID									
			}

			/**********************************************************************************************/
			// PREPARE FOR THE START OF THE MAIN LOOP
			// The proposed basis needs XtX and XtY of the current basis to compute XtX_prop and XtY_prop.
			// Other varialbes (e.g., beta and cholXtX) are not cross-used by the basis and basis_prop
			/**********************************************************************************************/
			U32 ite            = 0;
			U32 sample         = 0;

			/*************************************************************/
			// samplesInserted must be reset to zero at the start of each chain
			/*************************************************************/
			if (extra.computeCredible) { 
				for (int i = 0; i < NumCIVars; i++) {
					ci[i].samplesInserted = 0;
				} 
			} 

			PROP_DATA PROPINFO = { .N = N, .Npad16 = Npad16,  .samples = &sample, .keyresult = coreResults, .mem = Xnewterm,.model = &MODEL,
							  .pRND = &RND, .yInfo = &yInfo,  .sigFactor = opt->prior.sigFactor, .outlierSigFactor = opt->prior.outlierSigFactor,
							  .nSample_DeviationNeedUpdate = 1L, .shallUpdateExtremVec=0L, 
				               .numBasisWithoutOutlier=MODEL.NUMBASIS - (opt->prior.basisType[MODEL.NUMBASIS - 1] == OUTLIERID),};
	
			// Moved here bvz its xcols has two fixed lements, N and Npad
			NEWTERM   NEW = { .newcols = {.N=N, .Nlda=Npad} };
			
			I32 numBadIterations = 0;
			while (sample < MCMC_SAMPLES)
			{
				ite++;

				/**********************************************************************/
				/*     Re-generate a pool of random numbers if almost used up          */
				/***********************************************************************/
				// vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, MAX_RAND_NUM, rnd32, 0.f, 1.0f);
				if (RND.rnd32    >= RND32_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, (RND.rnd32 - RND32),                                 (U32PTR)RND32 ); RND.rnd32    = RND32;   }
				if (RND.rnd16    >= RND16_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, ((char*)RND.rnd16 - (char*)RND16 + 3) / sizeof(U32), (U32PTR)RND16 ); RND.rnd16    = RND16;   }
				if (RND.rnd08    >= RND08_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, ((char*)RND.rnd08 - (char*)RND08 + 3) / sizeof(U32), (U32PTR)RND08 ); RND.rnd08    = RND08;   }
				if (RND.rndgamma >= RNDGAMMA_END) {r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE,      stream, MAX_RAND_NUM, RNDGAMMA, (hyperPar.alpha_1+yInfo.n*0.5f), 0.f, 1.f  ); RND.rndgamma = RNDGAMMA;}
	
				// CHOOSE A BASI TYPE
				BEAST2_BASIS_PTR basis = MODEL.b + MODEL.PickBasisID(&PROPINFO); //basisID = *RND.rnd08++ < 128;				
				
				// IMPLEMENT THE NEW PROPOSED BASIS	
				// CHANGE: new.newKnot, numSeg, SEG/R1/R2, orders2, newIdx, nKnot_new, jumpType,propinfo.pRND.rnd8/32	
				// k1_old and k2_old in NEW.newcols are updated. And PROPINFO->mem=Xnewterm is used as a temp membuf here
				basis->Propose(basis, &NEW, &PROPINFO);

				#if  DEBUG_MODE ==1
					I32 basisIdx = basis - MODEL.b;		
					flagSat[NEW.newKnot - 1] += basisIdx == 0 && (NEW.jumpType == BIRTH || NEW.jumpType == MOVE);					 
					if (basisIdx == 0) ++(flagS[NEW.jumpType]);
					else 		       ++(flagT[NEW.jumpType]); 

					MEM.verify_header(&MEM);
				#endif

				/**********************************************************************/
				// Generate the new terms for the propsed step: To add or move a bk, two segments 
				// are re-generated; to remove a bk or merge two bks, one segment are re-generated.	
				/**********************************************************************/								
				I32 Knewterm = 0;				
				for (I32 i = 0; i < NEW.numSeg; i++) { 
					// NEW.numSeg is 0 if removing terms for ChORDER(newOrder<oldeTerm)
					// [1]..(1)....(2)...(sY-1)...N [N+1] 
					I32 kterms   = NEW.SEG[i].K = basis->GenTerms(Xnewterm+Npad*Knewterm, N, &NEW.SEG[i], &(basis->bConst));					
					Knewterm    += kterms;
				} // Iterate through all the new segments

				NEW.newcols.k2_new = NEW.newcols.k1 + Knewterm - 1L;	// if Knewterm=0 (i.e., delete terms), k2_new < k1_new
				
				// Get k1_old, k2_old, k1_new and k2_new by adding the base start
				NEW.newcols.k1     += basis->Kbase;		// k1=k1_new=k1_old
				NEW.newcols.k2_old += basis->Kbase;
				NEW.newcols.k2_new += basis->Kbase;
					 
				int KOLD = MODEL.curr.K;
				int KNEW = MODEL.curr.K + NEW.newcols.k2_new - NEW.newcols.k2_old;
				NEW.newcols.Knewterm = Knewterm;
				NEW.newcols.KOLD     = KOLD; // Total number of basis for the current model	
				NEW.newcols.KNEW     = KNEW; // Total number of bases in the proposed model	    

				/*************************************************************************/
				// Get XtX_prop: Copy parts of XtX to XtT_prop and fill new components 
				/*************************************************************************/
				// Set those rows of Xt_mars_newterms specfied by rowsMissing  to zeros
				if (yInfo.nMissing > 0 && Knewterm > 0 /*&& basis->type != OUTLIERID*/)  // needed for basisFunction_OUliter=1
				f32_mat_multirows_extract_set_by_scalar(Xnewterm,Npad,Knewterm,Xt_zeroBackup, yInfo.rowsMissing, yInfo.nMissing, 0.0f);

			    #if DEBUG_MODE == 1
					MEM.verify_header(&MEM);
				#endif

				/*********************************************************************************/
				//           Compute XtX_rop and XtY_prop from XtX andXtY
				/********************************************************************************/
	 			Update_XtX_from_Xnewterm(Xt_mars, Xnewterm, MODEL.curr.XtX, MODEL.prop.XtX, &NEW, &MODEL); // MODEL not used if GroupMat==1
				Update_XtY_from_Xnewterm(yInfo.Y, Xnewterm, MODEL.curr.XtY, MODEL.prop.XtY, &NEW, q);

				/*************************************************************/
				//  XtX_prop has been constructed, then compute chol_XtX_prop
				/*************************************************************/
				if (1L) {
					//Add precison values to the diagonal of XtX: post_P=XtX +diag(prec) 
					//Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
					/*				
					//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
					r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', KNEW, MODEL.prop.cholXtX, KNEW); // Choleskey decomposition; only the upper triagnle elements are used
					*/

					for (I32 i = 1; i < NEW.newcols.k1; i++) 	SCPY(i, MODEL.curr.cholXtX + (i - 1) * KOLD, MODEL.prop.cholXtX + (i - 1) * KNEW);

					precFunc.SetPropPrecXtXDiag_NtermsPerGrp(&MODEL, basis, &NEW); //&NEW is used only for OrderWise
					precFunc.chol_addCol(MODEL.prop.XtX + (NEW.newcols.k1 - 1) * KNEW, MODEL.prop.cholXtX, MODEL.prop.precXtXDiag, KNEW, NEW.newcols.k1, KNEW);

					/*
					for (rI32 i = 1; i <= (NEW.k1_new - 1L); i++) 	r_cblas_scopy(i, MODEL.curr.cholXtX + (i - 1L) * KOLD, 1L, MODEL.prop.cholXtX + (i - 1L) * KNEW, 1L);
					chol_addCol(MODEL.prop.cholXtX+ (NEW.k1_new - 1L)*KNEW, MODEL.prop.cholXtX, KNEW, NEW.k1_new, KNEW);
	
					{
					//LAPACKE_dpotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
					 SCPY(KNEW, MODEL.prop.XtY, MODEL.prop.beta_mean);
					 r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', KNEW, 1, MODEL.prop.cholXtX, KNEW, MODEL.prop.beta_mean, KNEW);
					}*/
				}
				
				/****************************************************************************************/
				/*    Compute Marg)like for the proposed model: prop.K is needed for computing prop.marg_Lik                                     */
				/****************************************************************************************/

				// In MODEL, basis's K is still the old ones. They will be updated only if the proposal is accetped
				// via basis->CalcBasisKsKeK_TermType(basis). The Ks of the proposed bases are stored in MODEL.prop.ntermP				   				/***************************************************************/	
				MODEL.prop.K = KNEW;
				precFunc.ComputeMargLik(&MODEL.prop, &MODEL.precState, &yInfo, &hyperPar);
				//if (MODEL.prop.marg_lik != MODEL.prop.marg_lik || fabs(MODEL.prop.marg_lik )>FLOAT_MAX || MODEL.prop.alpha2_star <0.f) {					   				  

				if (IsNaN(MODEL.prop.marg_lik) || IsInf(MODEL.prop.marg_lik)) {
					if (++numBadIterations < 15) {
						precFunc.IncreasePrecValues(&MODEL);
						precFunc.SetPrecXtXDiag(MODEL.curr.precXtXDiag, MODEL.b, MODEL.NUMBASIS, &MODEL.precState);
						precFunc.chol_addCol(MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.curr.precXtXDiag, MODEL.curr.K, 1L, MODEL.curr.K);
						precFunc.ComputeMargLik(&MODEL.curr, &MODEL.precState, &yInfo, &hyperPar);

						
					#if !( defined(R_RELEASE) || defined(M_RELEASE) || defined(P_RELEASE) )
						r_printf("prec: %.4f| marg_lik_prop: %.4f | marg_like_curr: %.4f \n", MODEL.precState.precVec[0], MODEL.prop.marg_lik, MODEL.curr.marg_lik);
					#endif
						continue;
					} else {
						skipCurrentPixel = 2;
						break;
					}

				}	else {
					numBadIterations = 0;
				} //if (marg_lik_prop != marg_lik_prop || alpha2_star_prop <0.f) 
 

			   /****************************************************************************************/
			   /*    DETERMINE WHETHER OR NOT TO ACCCEPT THE PROPSOED STEP                              */
			   /****************************************************************************************/
						
				F32 delta_lik = MODEL.prop.marg_lik - MODEL.curr.marg_lik;
				if   ( !(NEW.jumpType ==MOVE || basis->type ==OUTLIERID) ) {
					// Calcuate a factor adjusting the likelihood change
					F32 factor = basis->ModelPrior(basis, &NEW.newcols, &NEW); 
					delta_lik += factor;
				}
				
				//acceptTheProposal = *(RND.rnd16)++ < fastexp(delta_lik) * 65535.0f;				
				I08     acceptTheProposal;
				if      (delta_lik >   0.0f)   acceptTheProposal  = 1;
				else if (delta_lik < -23.0f)   acceptTheProposal  = 0;				
				else {				 
					F32    expValue = fastexp(delta_lik);
					if     (delta_lik > -0.5f)    acceptTheProposal = *(RND.rnd08)++ < expValue * 255.0f;
					else if(delta_lik > -5.0f  )  acceptTheProposal = *(RND.rnd16)++ < expValue * 65535.0f;
					else                          acceptTheProposal = *(RND.rnd32)++ < expValue * 4.294967296e+09;					 
				}
		 
				#if DEBUG_MODE == 1
                     if (basisIdx == 0) ++(flagS[NEW.jumpType]);
                     else 		   ++(flagT[NEW.jumpType]);
                     MEM.verify_header(&MEM);
				#endif

				if(acceptTheProposal)
				{
					#if DEBUG_MODE == 1
						if (basisIdx == 0) ++(accS[NEW.jumpType]);
						else               ++(accT[NEW.jumpType]);
					#endif

					//Recover the orignal vaules for those rows corresponding to missing Y values
					if (yInfo.nMissing > 0 && Knewterm > 0 /*&& basis->type != OUTLIERID*/)  //needed for basisFunction_OUliter=1						
						f32_mat_multirows_set_by_submat(Xnewterm, Npad, Knewterm, Xt_zeroBackup, yInfo.rowsMissing, yInfo.nMissing);

					// Inserting XnewTerms into Xt_mars
					if (NEW.newcols.k2_old != KOLD && NEW.newcols.k2_new != NEW.newcols.k2_old)
						shift_lastcols_within_matrix(Xt_mars, Npad, NEW.newcols.k2_old+1, KOLD, NEW.newcols.k2_new+1);
					if (Knewterm != 0)
						SCPY(Knewterm*Npad, Xnewterm, Xt_mars + (NEW.newcols.k1-1) * Npad);
					
					/****************************************************/
					//Find the good positions of the proposed MOVE, then update the knotLists and order
					/****************************************************/
 				 
					// The temporay changes made to KNOT are used for TREND and SEASON bases not OUTLIER: 
					 (basis->KNOT[-1] = basis->KNOT[INDEX_FakeStart], basis->KNOT[basis->nKnot]=basis->KNOT[INDEX_FakeEnd]);
					 basis->UpdateGoodVec_KnotList(basis, &NEW, Npad16);
					 (basis->KNOT[-1] = 1,   				          basis->KNOT[basis->nKnot] = N + 1L);

					basis->CalcBasisKsKeK_TermType(basis);
					UpdateBasisKbase(MODEL.b, MODEL.NUMBASIS, basis-MODEL.b);//basisIdx=basis-b Re-compute the K indices of the bases after the basisID 
	
					// Switching between Basis and Basis_prop
					{   //http: //stackoverflow.com/questions/3647331/how-to-swap-two-numbers-without-using-temp-variables-or-arithmetic-operations
						//basis  = ((I64)basis ^ (I64)basis_prop);basis_prop = ((I64)basis ^ (I64)basis_prop); basis      =  ((I64)basis ^ (I64)basis_prop);						
					
						#define Exchange(x)          {VOIDPTR tmp=MODEL.curr.x;  MODEL.curr.x=MODEL.prop.x; MODEL.prop.x=tmp;}
						#define Exchange2(x,y)       Exchange(x) Exchange(y) 
						#define Exchange4(x,y,z,k)   Exchange(x) Exchange(y) Exchange(z) Exchange(k)

						Exchange4(XtX, XtY, cholXtX, beta_mean);   // No need to exchange beta
						Exchange2(precXtXDiag, nTermsPerPrecGrp);  // needed for compontwise and orderwise						
						Exchange(alpha2Q_star);                    // changed for MRBEAST 						
						MODEL.curr.marg_lik    = MODEL.prop.marg_lik;
						MODEL.curr.K           = MODEL.prop.K;     //GetNumOfXmarCols(&MODEL): this function should also give KNEW; if not, there must be something wrong!						
					}

					#if  DEBUG_MODE ==1 
					if (q == 1) {
						//BEAST2_EvaluateModel(&MODEL.prop, MODEL.b, Xdebug, N, MODEL.NUMBASIS, &yInfo, &hyperPar, MODEL.precState.precVec, &stream);											
						//r_printf("ite:%d K: |%f|%f|diff:%f\n", ite, MODEL.curr.K, MODEL.curr.marg_lik, MODEL.prop.marg_lik, MODEL.prop.marg_lik - MODEL.curr.marg_lik);
					    //r_printf(" %f[%f]-%f %f\n", (MODEL.curr.alpha2_star), (MODEL.prop.alpha2_star),	yInfo.alpha1_star* (log(MODEL.curr.alpha2_star) - log(MODEL.prop.alpha2_star)), yInfo.alpha1_star);
						
						/*
						I32 K = MODEL.prop.K;
						for (int i = 0; i < K; ++i) {

							r_printf("%f %f %f | %f %f %f\n",
								MODEL.prop.cholXtX[i * K + i], MODEL.curr.cholXtX[i * K + i],
								MODEL.prop.cholXtX[i * K + i] - MODEL.curr.cholXtX[i * K + i],
								MODEL.prop.XtX[i * K + i], MODEL.curr.XtX[i * K + i],
								MODEL.prop.XtX[i * K + i] - MODEL.curr.XtX[i * K + i]);

						}
						r_printf("ite----%d\n", ite);
						int a = 1;
						 */
					}
					else {
						BEAST2_EvaluateModel(&MODEL.prop, MODEL.b, Xdebug, N, MODEL.NUMBASIS, &yInfo, &hyperPar, MODEL.precState.precVec, &stream);
						//r_printf("MRite%d |%f|%f|diff:%f -prec %f\n", ite, MODEL.curr.marg_lik, MODEL.prop.marg_lik, MODEL.prop.marg_lik - MODEL.curr.marg_lik, MODEL.precState.precVec[0]);
					 
	                                   /****
						I32 K = MODEL.prop.K;
						for (int i = 0; i < MODEL.prop.K; ++i) {
						 
								r_printf("%f %f %f | %f %f %f\n", 
									MODEL.prop.cholXtX[i*K + i], MODEL.curr.cholXtX[i * K + i],
									MODEL.prop.cholXtX[i * K + i] - MODEL.curr.cholXtX[i * K + i],
									MODEL.prop.XtX[i * K + i], MODEL.curr.XtX[i * K + i],
									MODEL.prop.XtX[i * K + i] - MODEL.curr.XtX[i * K + i]);
						 
						}
						
							r_printf("ite----%d\n",ite);
							int a = 1;
					   ****/ 
					}

					#endif

				} //(*rnd32++ < exp(marg_lik_prop - basis->marg_lik))
				
				/****************************************************************************************/
				// For buin-in iterations, we also want to sample parameters: added for cases when burn-in is very large
				/****************************************************************************************/
				const int bResampleParameter  = (ite >=100)         && (ite % 20 == 0)     ;      //MCMC_BURIN>=150L
				const int bStoreCurrentSample = (ite > MCMC_BURNIN) && (ite % MCMC_THINNING == 0);

				 /*******************************************/
				// Re-SAMPLING SIG2,  Beta, and PrecValue
				/*******************************************/

				if (q == 1) {
			 
				    //First, Re-SAMPLING SIG2			 
					if (bResampleParameter || bStoreCurrentSample)	{	 
						//vdRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, &sig2, (alpha_1+n/2), 0, 1.0/(alpha_1+basis->alpha2_star *0.5));
						F32 sig2_inv  = (*RND.rndgamma++) * 1.f / MODEL.curr.alpha2Q_star[0];
						F32 sig2      = 1.0f / sig2_inv;
						MODEL.sig2[0]     = sig2 > MIN_SIG2_VALUE ? sig2 : MODEL.sig2[0];
						//r_printf("ite-%d SIG %f %f\n", ite, MODEL.sig2,  MODEL.sig2*yInfo.sd*yInfo.sd);
					}	 

					//Re-sample beta to be used for either re-sampling prec (ite%20=0) or predicting Y (ite%thiningFactor=0)
					if (bResampleParameter || (bStoreCurrentSample && extra.useRndBeta)) {
							//Compute beta = beta_mean + Rsig2 * randn(p, 1);
							//Usig2 = (1 / sqrt(sig2)) * U; 		beta = beta_mean + linsolve(Usig2, randn(p, 1), opts);
							//status = vdRngGaussian( method, stream, n, r, a, sigma );
							I32 K = MODEL.curr.K;
							r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K, MODEL.beta, 0, 1);
							//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, MODEL.curr.cholXtX, K, MODEL.curr.beta, K); // LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
							solve_U_as_U_invdiag(MODEL.curr.cholXtX,  MODEL.beta, K, K);
							r_ippsMulC_32f_I(fastsqrt(MODEL.sig2[0]), MODEL.beta, K);
							r_ippsAdd_32f_I(MODEL.curr.beta_mean,     MODEL.beta, K);
					}

					// Re-sample the precison parameters and re-calcuate marg_lik and beta
					if (bResampleParameter ) {
						/*
						 I32   K     = MODEL.K;
						 F32   sumq  = DOT(K, MODEL.curr.beta, MODEL.curr.beta);
						 r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec, (hyperPar.del_1 + K * 0.5f), 0, 1.f);
						 modelPar.prec[2]     = modelPar.prec[1] = modelPar.prec[0] = (*modelPar.prec) / (hyperPar.del_2 + 0.5f * sumq / MODEL.sig2);
						 modelPar.LOG_PREC[2] = modelPar.LOG_PREC[1] = modelPar.LOG_PREC[0] = logf(modelPar.prec[0]);
						 */

						 /*
						 //X_mars_prop has been constructed. Now use it to calcuate its marginal likelihood
						 //Add precison values to the diagonal of XtX: post_P=XtX +diag(prec)
						 SCPY(K * K, MODEL.curr.XtX, MODEL.curr.cholXtX);

						 {//Add precison values to the SEASONAL diagonal compoents
							 //U08PTR termType = MODEL.termType;
							 for (I32 i = 1, j = 0; i <= K; i++)						{
								 //MODEL.curr.cholXtX[j + (i)-1] += modelPar.prec[*termType++];
								 MODEL.curr.cholXtX[j + (i)-1] += modelPar.prec[0];
								 j += K;
							 }
						 }//Add precison values to the SEASONAL diagonal compoents

						 //Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
						 //lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
						 r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', K, MODEL.curr.cholXtX, K); // Choleskey decomposition; only the upper triagnle elements are used
						 //chol_addCol(MODEL.curr.cholXtX, MODEL.curr.cholXtX, K, 1, K);

						 //LAPACKE_spotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
						 SCPY(K, MODEL.curr.XtY, MODEL.curr.beta_mean);
						 r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', K, 1, MODEL.curr.cholXtX, K, MODEL.curr.beta_mean, K);
						 */

						I32 ntries = 0;
						do {
							if (ntries++ == 0)	precFunc.ResamplePrecValues(&MODEL, &hyperPar, &stream);
							else				precFunc.IncreasePrecValues(&MODEL);
							precFunc.SetPrecXtXDiag(MODEL.curr.precXtXDiag, MODEL.b, MODEL.NUMBASIS, &MODEL.precState);
							precFunc.chol_addCol(MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.curr.precXtXDiag, MODEL.curr.K, 1L, MODEL.curr.K);
							precFunc.ComputeMargLik(&MODEL.curr, &MODEL.precState, &yInfo, &hyperPar);
						} while (IsNaN(MODEL.curr.marg_lik) && ntries < 20);

						if (IsNaN(MODEL.curr.marg_lik)) {
                             #if !(defined(R_RELEASE) || defined(M_RELEASE) ||  defined(P_RELEASE)) 
							 r_printf("skip3 | prec: %.4f| marg_lik_cur: %.4f \n", MODEL.precState.precVec[0], MODEL.curr.marg_lik);
                             #endif
							skipCurrentPixel = 3;
							break;
						}

						/* No need to re-sample beta because it is not really used
						r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K, beta, 0, 1);
						// LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
						r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, cholXtX, K, beta, K);
						r_ippsMulC_32f_I(sqrtf(modelPar.sig2), beta, K);
						r_ippsAdd_32f_I(beta_mean, beta, K);
						*/
						/* /FInally, re-sample sig2 based on the lastest alpha2_star
						//vdRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, &sig2, (alpha_2+n/2), 0, 1.0/(alpha_1+basis->alpha2_star *0.5));
						modelPar.sig2 = (*rndgamma++)*1.0f / (modelPar.alpha_1 + basis->alpha2_star *0.5f);
						modelPar.sig2 = 1.f / modelPar.sig2;
						*/
					}

				} // if (q == 1)

				// For MRBEAST only
				if (q != 1) {    

					F32PTR MEMBUF = Xnewterm;
					I32    K      = MODEL.curr.K;
					//First, Re-SAMPLING SIG2	 
					if (bResampleParameter || bStoreCurrentSample) {						
						local_pcg_invwishart_upper(&stream, MODEL.sig2, MODEL.sig2 + q * q, MEMBUF, q, MODEL.curr.alpha2Q_star, hyperPar.alpha_1 + yInfo.n + q - 1);
					}

					// RE-SAMPLE beta to be used for either re-sampling prec (ite%20=0) or predicting Y (ite%thiningFactor=0)
					if (bResampleParameter || (bStoreCurrentSample && extra.useRndBeta)) {
						r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K * q, MEMBUF, 0., 1.);
						r_cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, K, q, q, 1.0, MEMBUF, K, MODEL.sig2, q, 0.f, MODEL.beta, K);
						//LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
						//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, basis->post_P_U, K, beta, K);
						//r_ippsMulC_32f_I(sqrtf(modelPar.sig2), beta, K);
						//r_ippsAdd_32f_I(basis->beta_mean, beta, K);

						//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, q, post_P_U, K, beta, K);
						//r_ippsAdd_32f_I(MODEL.curr.beta_mean, MODEL.beta, K * q);
						solve_U_as_U_invdiag_multicols(MODEL.curr.cholXtX, MODEL.beta, K, K, q);
						r_ippsAdd_32f_I(MODEL.curr.beta_mean, MODEL.beta, K * q);
					}

					// Re-sample the precison parameters and re-calcuate marg_lik and beta
					if (bResampleParameter )		{						
						//FLOAT_SHARE.sumq = DOT(K, basis->beta, basis->beta);
						//r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec, (modelPar.alpha_2 + K *0.5f), 0, 1.f);
						//modelPar.prec[2] = modelPar.prec[1] = (*modelPar.prec) / (modelPar.alpha_1 + 0.5f*FLOAT_SHARE.sumq / modelPar.sig2);
						I32 ntries = 0;
						do {
							if (ntries++ == 0) {
								// Get trace( B*inv(SIG2)*B')						
								//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', q, q, basis->alpha_Q_star, q, W_L, q);	
								r_cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, K, q, q, 1.0, MODEL.beta, K, MODEL.sig2 + q * q, q, 0.f, MEMBUF, K);
								F32 sumq = DOT(K * q, MEMBUF, MEMBUF);
								r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1L, MODEL.precState.precVec, (hyperPar.del_1 + K * q * 0.5f), 0.f, 1.f);
								MODEL.precState.precVec[0]    = MODEL.precState.precVec[0] / (hyperPar.del_2 + 0.5f * sumq);
								MODEL.precState.logPrecVec[0] = logf(MODEL.precState.precVec[0]);
							}	else {
								precFunc.IncreasePrecValues(&MODEL);
							}
							//precFunc.ResamplePrecValues( &MODEL, &hyperPar,&stream);
							precFunc.SetPrecXtXDiag(MODEL.curr.precXtXDiag, MODEL.b, MODEL.NUMBASIS, &MODEL.precState);
							precFunc.chol_addCol(MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.curr.precXtXDiag, K, 1, K);
							precFunc.ComputeMargLik(&MODEL.curr, &MODEL.precState, &yInfo, &hyperPar);
						} while (IsNaN(MODEL.curr.marg_lik) && ntries < 20);


						if (IsNaN(MODEL.curr.marg_lik)) {
                            #if !(defined(R_RELEASE) || defined(M_RELEASE)  || defined(P_RELEASE))
							 r_printf("skip4 | prec: %.4f| marg_lik_cur: %.4f \n", MODEL.precState.precVec[0], MODEL.curr.marg_lik);
                            #endif
							skipCurrentPixel = 3;
							break;
						}

					} // if (bResampleParameter )


				} // 	if (q != 1) {    

  	 
 
				if (!bStoreCurrentSample) continue;
				
				#if DEBUG_MODE == 1
					MEM.verify_header(&MEM);
				#endif

				/**********************************************/
				//
				//      Start to compute final results
				//
				/**********************************************/

				sample++;

				*resultChain.marg_lik += MODEL.curr.marg_lik;

				if (q == 1) {
					resultChain.sig2[0] += MODEL.sig2[0];
				}	else {
					// For MRBEAST
					F32PTR MEMBUF = Xnewterm;
					r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, q, 1.f, MODEL.sig2, q, MODEL.sig2, q, 0.f, MEMBUF, q);
					r_ippsAdd_32f_I(MEMBUF, resultChain.sig2,  q*q);
				}


				const F32PTR BETA = (extra.useRndBeta == 0) ? MODEL.curr.beta_mean : MODEL.beta;
				{
					 F32PTR MEMBUF1 = Xnewterm;

					for (I32 i = 0; i < MODEL.NUMBASIS; ++i) 
					{
						const BEAST2_BASIS_PTR  basis   = MODEL.b + i;
						const CORESULT        * result  = coreResults+i;
						const I32        nKnot  = basis->nKnot;
						const TKNOT_PTR  KNOT   = basis->KNOT;

						result->xNProb[nKnot] += 1L;

						//Counting probability of being breakpoints				
						for (I32 j = 0; j < nKnot; j++) result->xProb[ KNOT[j]-1 ] += 1L;

						//Summng up the harmonic orders or trend orders for individual seeasonal segments
						if (result->xorder != NULL) {
							for (I32 j = 0; j <= nKnot; ++j) {
								I32 r1 = KNOT[j-1], r2 = KNOT[j]-1;
								r_ippsAddC_32s_ISfs(basis->ORDER[j], result->xorder+r1 - 1, r2 - r1 + 1, 0);
							}
						}

						//Compute the averaged  signals
						if (q == 1) {
							//r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K_SN, 1.f, Xt_mars, Npad, MODEL.curr.beta_mean, 1L, 0.f, MEMBUF1, 1L);
							basis->ComputeY(Xt_mars, BETA, MEMBUF1, basis, Npad);
							f32_add_v_v2_vec_inplace(MEMBUF1, result->x, result->xSD, N);
							MEMBUF1 += Npad;
						}	else {
							// for MRBEAST
							F32PTR 	X    = Xt_mars + basis->Kbase * Npad;
							F32PTR  beta = BETA+basis->Kbase;
							I32     K    = basis->K;
							//r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, X, Npad, beta, 1L, 0.f, Y, 1L);
							r_cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, q, K, 1.0f,
                                                                       X, Npad, beta,  MODEL.curr.K, 0.f, MEMBUF1, N);
							f32_add_v_v2_vec_inplace(MEMBUF1, result->x, result->xSD, N*q);
							MEMBUF1 += Npad*q;
						}
	
						
					}//for (rI32 i = 0; i < MODEL.NUMBASIS; i++)
				
				}
				/*********************************************************************/
				// At this point, MEMBUF1=season, MEMBUF1+Npad=trend, and MEMBUF1+Npad+Npad=outlier;
				/*********************************************************************/

				/********************************************/
				// Compute results for the seasonal cmpnt
				/********************************************/
				if(extra.computeSeasonAmp) 
				{
					const F32PTR           MEMBUF1 = Xnewterm + 3*Npad;
	
					const BEAST2_BASIS_PTR basis    = &MODEL.b[MODEL.sid];
					const I32              knotNum  = basis->nKnot;
					const TKNOT_PTR        knotList = basis->KNOT;
					
					//Summng up the per-segment harmonic magnitudes  	
					F32PTR             beta            = BETA;
					const TORDER_PTR   orderList       = basis->ORDER;
					const F32PTR       SEASON_SQR_CSUM = basis->bConst.season.SQR_CSUM +1L;  //SQR_CSUM has a row length of (N+1)
					const F32PTR       SEASON_SCALE    = basis->bConst.season.SCALE_FACTOR;

					if (SEASON_BTYPE == 0) {
						for (I32 i = 0; i <= knotNum; i++) {
							const I32	r1       = knotList[i - 1];
							const I32	r2       = knotList[i] - 1;							
							const I32   order2   = orderList[i] * 2L;
							F32PTR seasonSqrCsum = SEASON_SQR_CSUM;							
							F32    amp           = 0;
							for (I32 j = 0; j < order2; j++) {
								//TODO: re-check here
								F32 scalingFactor = N / (seasonSqrCsum[r2 - 1] - seasonSqrCsum[(r1 - 1) - 1]);
								F32 beta0         = beta[j]* SEASON_SCALE[j];
								amp               = amp + (beta0 * beta0) * scalingFactor;
								seasonSqrCsum     += (N + 1LL);
							}			 
							//r_ippsSubC_32f_I(-amp, resultChain.samp + r1 - 1, segLength, 0); 
							r_ippsSet_32f(amp, MEMBUF1 + r1 - 1, r2 - r1 + 1L);
							beta += order2;
						}
					} else {
						for (I32 i = 0; i <= knotNum; i++) {
							const I32	r1       = knotList[i - 1];
							const I32	r2       = knotList[i] - 1;							
							const I32   order2   = orderList[i] * 2L;			
							F32    amp           = 0;		 
							for (I32 j = 0; j < order2; ++j) {
								F32    beta0 = beta[j] * SEASON_SCALE[j];
								amp += beta0 * beta0;
							}					 

							//r_ippsSubC_32f_I(-amp, resultChain.samp + r1 - 1, segLength, 0); 
							r_ippsSet_32f(amp, MEMBUF1 + r1 - 1, r2 - r1 + 1L);
							beta += order2;
						}
					}
	
					r_ippsAdd_32f_I(MEMBUF1, resultChain.samp,   N);
					r_ippsMul_32f_I(MEMBUF1, MEMBUF1,            N);
					r_ippsAdd_32f_I(MEMBUF1, resultChain.sampSD, N); //added to the square of the samp for computering SD

					
					if (extra.tallyPosNegSeasonJump) { 			 
						I32  posKnotNum = 0;
						for (I32 i = 0; i < knotNum; i++) { // It must be i<KnotNum bcz of dealing with knots only 
							I64 knot = knotList[i];
							if (MEMBUF1[knot - 1] > MEMBUF1[knot - 1 - 1]) {
								resultChain.spos_cpOccPr[knot - 1] += 1;
								posKnotNum++;
							}
						}
						resultChain.spos_ncpPr[posKnotNum]           += 1L;
						resultChain.sneg_ncpPr[knotNum - posKnotNum] += 1L;
					}
				}
				
				/********************************************/
				// Compute results for the trend cmpnt
				/********************************************/
				if(extra.computeTrendSlope)
				{
					const BEAST2_BASIS_PTR basis    = &MODEL.b[MODEL.tid];
					const I32              knotNum  = basis->nKnot;
					const TKNOT_PTR        knotList = basis->KNOT;

					const F32PTR TREND = Xnewterm + Npad * MODEL.tid;      //trend signal, already filled with real values
					const F32PTR SLP   = Xnewterm + Npad * MODEL.NUMBASIS; //slope: to be computed

																	// Compute the rate of change in trend based on beta. 
					f32_diff_back(TREND, SLP, N);
					f32_add_v_v2_vec_inplace(SLP, resultChain.tslp, resultChain.tslpSD, N); //added to the square of the trend signal for computing SD
					//i32_increment_bycond_inplace(resultChain.tslpSgnPosPr, SLP, N); //increamnent tslpSingPr if SLP is larger than 0				 
					i32_increment_vec2_bycond_inplace(resultChain.tslpSgnPosPr, resultChain.tslpSgnZeroPr, SLP, N); //increamnent tslpSingPr if SLP is larger than 0				 
					if (extra.tallyPosNegTrendJump ) {//NEWLY ADDEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD					 
						I32  posKnotNum = 0;
						for (I32 i = 0; i < knotNum; i++) {  // It must be i<sKnotNum bcz of dealing with knots only 
							I64 knot = knotList[i];
							if (SLP[knot - 1] > 0) {
								resultChain.tpos_cpOccPr[knot - 1] += 1;
								posKnotNum++;
							}
						}
						resultChain.tpos_ncpPr[posKnotNum]           += 1L;
						resultChain.tneg_ncpPr[knotNum - posKnotNum] += 1L;
					}

					if (extra.tallyIncDecTrendJump ){			 						
						I32  incKnotNum = 0;
						for (I32 i = 0; i < knotNum; i++) {  // It must be i<sKnotNum bcz of dealing with knots only 
							I64 knot = knotList[i];
							if (knot >= 2 && SLP[(knot + 1) - 1] > SLP[(knot - 1) - 1]) {
								resultChain.tinc_cpOccPr[knot - 1] += 1;
								incKnotNum++;
							}
						}
						resultChain.tinc_ncpPr[incKnotNum]         += 1L;
						resultChain.tdec_ncpPr[knotNum-incKnotNum] += 1L;
					}
			
				}

				/********************************************/
				// Compute results for the outlier cmpnt
				/********************************************/
				if(extra.tallyPosNegOutliers) {

					const BEAST2_BASIS_PTR basis    = &MODEL.b[MODEL.oid];
					const I32              knotNum  = basis->nKnot;
					const TKNOT_PTR        knotList = basis->KNOT;

					const F32PTR OUTLIIER  = Xnewterm + Npad* MODEL.oid;
	 
					I32  posKnotNum = 0;
					for (I32 i = 0; i < knotNum; i++) {  // It must be i<sKnotNum bcz of dealing with knots only 
						I64 knot = knotList[i];
						if (OUTLIIER[knot - 1] > 0) {
							resultChain.opos_cpOccPr[knot - 1] += 1;
							posKnotNum++;
						}							
					}
					resultChain.opos_ncpPr[posKnotNum]			 += 1L;
					resultChain.oneg_ncpPr[knotNum - posKnotNum] += 1L;								
				}

				/*************************************************/
				// Compute ci intervals: new row of data have already calculated
				// and saved in Xnewterm, Xnewterm+Npad, and Xnweterm+2*Npad
				/*************************************************/
				if (extra.computeCredible)	{ 	

					// when  *RND.rnd16++ <= ciParam.subsampleFraction_x_INT16MAX, the current sample not included.
					// otherwise, no need to insert it into the ci strips. So, just skip to the next iteration					
					if ( !extra.fastCIComputation ||   *RND.rnd16++  < ciParam.subsampleFraction_x_INT16MAX  ) {
						//if (*rnd32++ < subsampleFraction*4.294967296000000e+09)
	    
						// New row of data for  seasonal, trend, outelier, and slope
						for (int i = 0; i <  NumCIVars; i++) 
							InsertNewRowToUpdateCI(&ciParam, &ci[i]);						
					}					

				} // if (extra.computeCredible)
				

				// dumpMCMCsamples==1 only if q=1 and NumPiexl=1
				if (extra.dumpMCMCSamples) {
				
					F32PTR cmpnts = Xnewterm;
					
					if (result.smcmc) {
						f32_copy(cmpnts, result.smcmc, N*q);
						cmpnts       += N*q;
						result.smcmc += N * q;
					}

					f32_copy(cmpnts, result.tmcmc, N* q);
					cmpnts += N * q;
					result.tmcmc += N * q;

					if (result.omcmc) {
						f32_copy(cmpnts, result.omcmc, N * q);
						result.omcmc += N * q;
					}
				}

				if (extra.printProgress && NUM_PIXELS == 1 && sample % 1000 == 0) {
					F32 frac = (F32)(chainNumber * MCMC_SAMPLES + sample) / (MCMC_SAMPLES * MCMC_CHAINNUM);
					printProgress1(frac, extra.consoleWidth, Xnewterm, 0);
				}

			}//WHILE(sample<SAMPLE)

			/*****************************************************************************/
			//
			// One chain is done; then start post-processing the result of the chain     
			//
			/*****************************************************************************/

			if (!skipCurrentPixel)
			{
				int   sum;
				#define GetSum(arr) (r_ippsSum_32s_Sfs(arr, N, &sum, 0), sum) // parenthesis operator

				I32   sMAXNUMKNOT = MODEL.sid >=0 || MODEL.vid >= 0 ? MODEL.b[0].prior.maxKnotNum:-9999999; // the seasonal cmpt is always the first one
				I32   tMAXNUMKNOT = MODEL.tid>=0?  MODEL.b[MODEL.tid].prior.maxKnotNum:-9999999;
				I32   oMAXNUMKNOT = MODEL.oid>=0?  MODEL.b[MODEL.oid].prior.maxKnotNum:- 9999999;
				F32   inv_sample  = 1.f / sample;			
				

				*resultChain.marg_lik = *resultChain.marg_lik * inv_sample;
				// FOR MRBEAST
				for (int col = 0; col < q; col++)	{ 
					for (int i = 0; i < q; i++) {
						resultChain.sig2[col*q+i] = resultChain.sig2[col*q+i] * inv_sample * yInfo.sd[col] * yInfo.sd[i];
					}
				}
				
				///////////////SEASON////////////////////////////////////////////////////
				if (MODEL.sid >= 0 || MODEL.vid >= 0) {

						*resultChain.sncp = GetSum(resultChain.scpOccPr)* inv_sample; 
						i32_to_f32_scaleby_inplace(resultChain.sncpPr,	(sMAXNUMKNOT + 1), inv_sample);
						i32_to_f32_scaleby_inplace(resultChain.scpOccPr, N,		           inv_sample);	
					
						//FOR MRBEAST
						for (int i = 0; i < q; i++) {
							F32 offset = 0.0f;
							f32_sx_sxx_to_avgstd_inplace(resultChain.sY + i * N, resultChain.sSD + i * N, sample, yInfo.sd[i], offset, N);
						}

						if (extra.computeSeasonOrder) i32_to_f32_scaleby_inplace(resultChain.sorder, N, inv_sample);
						if (extra.computeSeasonAmp) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {
								F32 offset = 0.0f;
								f32_sx_sxx_to_avgstd_inplace(resultChain.samp+i*N, resultChain.sampSD + i*N, sample, yInfo.sd[i]* yInfo.sd[i], offset, N);
							}							
						}
						if (extra.computeCredible) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {								
								r_ippsMulC_32f_I(yInfo.sd[i], resultChain.sCI+      N*i, N);
								r_ippsMulC_32f_I(yInfo.sd[i], resultChain.sCI+N*q + N*i, N);
							}							
						} 
				}

				///////////////TREND////////////////////////////////////////////////////
				if (MODEL.tid >= 0) {
						*resultChain.tncp = GetSum(resultChain.tcpOccPr)*inv_sample ;
						i32_to_f32_scaleby_inplace(resultChain.tncpPr, (tMAXNUMKNOT + 1),	inv_sample);
						i32_to_f32_scaleby_inplace(resultChain.tcpOccPr, N,					inv_sample);
						//FOR MRBEAST
						for (int i = 0; i < q; i++) {
							F32 offset = yInfo.mean[i];
							f32_sx_sxx_to_avgstd_inplace(resultChain.tY + i * N, resultChain.tSD + i * N, sample, yInfo.sd[i],offset, N);
						}


						if (extra.computeTrendOrder) 	i32_to_f32_scaleby_inplace(resultChain.torder, N, inv_sample);	

						if (extra.computeTrendSlope) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {
								f32_sx_sxx_to_avgstd_inplace(resultChain.tslp+i*N, resultChain.tslpSD+ i*N, sample, yInfo.sd[i]/opt->io.meta.deltaTime, 0, N);
							}							
							i32_to_f32_scaleby_inplace(resultChain.tslpSgnPosPr, N*q, inv_sample);
							i32_to_f32_scaleby_inplace(resultChain.tslpSgnZeroPr, N*q, inv_sample);
						}

						if (extra.computeCredible) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {
								f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], resultChain.tCI +        N * i, N);
								f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], resultChain.tCI + N*q +  N * i, N);
								//r_ippsMulC_32f_I(,    resultChain.tCI+(2*N)*i, N + N),
								//r_ippsSubC_32f_I(-yInfo.mean[i], ); //ippsAddC_32f_I(yInfo.mean, result.tCI, N + N);
							}	

							if (extra.computeTrendSlope) {
								for (int i = 0; i < q; i++) {
									f32_mul_val_inplace(yInfo.sd[i],  resultChain.tslpCI + N * i, N);
									f32_mul_val_inplace(yInfo.sd[i],  resultChain.tslpCI + N * q + N * i, N);
									//r_ippsMulC_32f_I(,    resultChain.tCI+(2*N)*i, N + N),
									//r_ippsSubC_32f_I(-yInfo.mean[i], ); //ippsAddC_32f_I(yInfo.mean, result.tCI, N + N);
								}						
							}							
						}
						
				}

				///////////////OUTLIR////////////////////////////////////////////////////
				if (MODEL.oid >= 0) {					
					 *resultChain.oncp = inv_sample * GetSum(resultChain.ocpOccPr);
					 i32_to_f32_scaleby_inplace(resultChain.oncpPr,  (oMAXNUMKNOT + 1), inv_sample);
					 i32_to_f32_scaleby_inplace(resultChain.ocpOccPr, N,                inv_sample);
					 //FOR MRBEAST
					 for (int i = 0; i < q; i++) {
						 f32_sx_sxx_to_avgstd_inplace(resultChain.oY+i*N, resultChain.oSD+i*N, sample, yInfo.sd[i], 0, N);
					 }					 
					 if (extra.computeCredible) {
						 //FOR MRBEAST
						 for (int i = 0; i < q; i++) {
							 r_ippsMulC_32f_I(yInfo.sd[i], resultChain.oCI +      i*N, N );
							 r_ippsMulC_32f_I(yInfo.sd[i], resultChain.oCI + N*q+ i*N, N);
						 }						 
					 }	
				}


				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyPosNegSeasonJump && MODEL.sid>=0) {

					GetSum(resultChain.spos_cpOccPr);	
					*resultChain.spos_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.sneg_ncp = *resultChain.sncp - *resultChain.spos_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.spos_ncpPr, (sMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.sneg_ncpPr, (sMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.spos_cpOccPr, N, inv_sample);
 
					SCPY(N, resultChain.scpOccPr, resultChain.sneg_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.spos_cpOccPr, (F32PTR)resultChain.sneg_cpOccPr, N);   //NEWLY ADDED
				}

				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyPosNegTrendJump) {

					GetSum(resultChain.tpos_cpOccPr);	
					*resultChain.tpos_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.tneg_ncp = *resultChain.tncp - *resultChain.tpos_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.tpos_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tneg_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tpos_cpOccPr, N,				  inv_sample);
 
					SCPY(N, resultChain.tcpOccPr, resultChain.tneg_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.tpos_cpOccPr, (F32PTR)resultChain.tneg_cpOccPr, N);   //NEWLY ADDED
				}


				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyIncDecTrendJump) {

					GetSum(resultChain.tinc_cpOccPr);	   
					*resultChain.tinc_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.tdec_ncp = *resultChain.tncp - *resultChain.tinc_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.tinc_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tdec_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tinc_cpOccPr, N, inv_sample);

					SCPY(N, resultChain.tcpOccPr, resultChain.tdec_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.tinc_cpOccPr, (F32PTR)resultChain.tdec_cpOccPr, N);   //NEWLY ADDED
				}


				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyPosNegOutliers && MODEL.oid >= 0) {

					GetSum(resultChain.opos_cpOccPr);	
					*resultChain.opos_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.oneg_ncp = *resultChain.oncp - *resultChain.opos_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.opos_ncpPr, (oMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.oneg_ncpPr, (oMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.opos_cpOccPr, N, inv_sample);

					SCPY(N, resultChain.ocpOccPr, resultChain.oneg_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.opos_cpOccPr, (F32PTR)resultChain.oneg_cpOccPr, N);   //NEWLY ADDED
				}

			}// Finish computing the result of the single chain


			/**************************************************/
			//   Add up the individual chain to the Result
			/**************************************************/
			if (!skipCurrentPixel) 
			{
				I32   sMAXNUMKNOT = MODEL.sid >= 0 || MODEL.vid >= 0 ? MODEL.b[0].prior.maxKnotNum : -9999999;
				I32   tMAXNUMKNOT = MODEL.tid>=0?  MODEL.b[MODEL.tid].prior.maxKnotNum:-9999999;
				I32   oMAXNUMKNOT = MODEL.oid>=0?  MODEL.b[MODEL.oid].prior.maxKnotNum:- 9999999;

				//https://stackoverflow.com/questions/13216423/error-pasting-and-red-does-not-give-a-valid-preprocessing-token
				//#define _1(x)      *(result.##x) += *(resultChain.##x) // Working with MSVC but not GCC

				#define _1(x)      *(result.x) += *(resultChain.x)
				#define _N(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N)
				#define _Nq(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N*q)
				#define _q(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, q)
				#define _q2(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, q*q)
				#define _2N(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N+N)
				#define _2Nq(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N*q+N*q)
				#define _skn_1(x)  r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, sMAXNUMKNOT + 1)
				#define _tkn_1(x)  r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, tMAXNUMKNOT + 1)
				#define _okn_1(x)  r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, oMAXNUMKNOT + 1)

				_1(marg_lik);
				_q2(sig2);  //Fpr MRBEAST
				
				if (MODEL.sid >= 0 || MODEL.vid >0) {
					_1(sncp); _skn_1(sncpPr);	     _N(scpOccPr); _Nq(sY); _Nq(sSD);
					if (extra.computeSeasonOrder)    _N(sorder);
					if (extra.computeSeasonAmp)      _N(samp), _N(sampSD);
					if (extra.computeCredible)       _2Nq(sCI);
				}

				if (MODEL.tid >= 0) {					
					_1(tncp); _tkn_1(tncpPr);	   _N(tcpOccPr); _Nq(tY); _Nq(tSD);
					if (extra.computeTrendOrder)   _N(torder);
					if (extra.computeTrendSlope)   _N(tslp), _N(tslpSD),_N(tslpSgnPosPr), _N(tslpSgnZeroPr);
					if (extra.computeCredible)     _2Nq(tCI);
					if (extra.computeCredible && extra.computeTrendSlope ) _2Nq(tslpCI);
				}

				if (MODEL.oid >= 0) {
					_1(oncp); _okn_1(oncpPr);	_N(ocpOccPr); _Nq(oY); _Nq(oSD);
					if (extra.computeCredible)   _2Nq(oCI);
				}

				if (extra.tallyPosNegSeasonJump && MODEL.sid >=0) {
					_1(spos_ncp);         _1(sneg_ncp); //NEWLY ADDED					
					_skn_1(spos_ncpPr);   _skn_1(sneg_ncpPr); //NEWLY ADDED
					_N(spos_cpOccPr);     _N(sneg_cpOccPr); //NEWLY ADDED	
				}

				if (extra.tallyPosNegTrendJump ) {
					_1(tpos_ncp);            _1(tneg_ncp); //NEWLY ADDED					
				    _tkn_1(tpos_ncpPr); _tkn_1(tneg_ncpPr); //NEWLY ADDED
					_N(tpos_cpOccPr);         _N(tneg_cpOccPr); //NEWLY ADDED 
				}
				
				if (extra.tallyIncDecTrendJump) {
					_1(tinc_ncp);         _1(tdec_ncp); //NEWLY ADDED					
					_tkn_1(tinc_ncpPr); _tkn_1(tdec_ncpPr); //NEWLY ADDED
					_N(tinc_cpOccPr);     _N(tdec_cpOccPr); //NEWLY ADDED
				}
				
				if (extra.tallyPosNegOutliers && MODEL.oid >= 0) {
					_1(opos_ncp);            _1(oneg_ncp); //NEWLY ADDED					
					_okn_1(opos_ncpPr);      _okn_1(oneg_ncpPr); //NEWLY ADDED
					_N(opos_cpOccPr);        _N(oneg_cpOccPr); //NEWLY ADDED 
				}	
 
				#undef _1
				#undef _N
				#undef _Nq 
				#undef _q 
				#undef _q2
				#undef _2N
				#undef _2Nq
				#undef _skn_1
				#undef _tkn_1
			    #undef _okn_1
			}


			// Jump out of the chainumber loop
			if (skipCurrentPixel) {
			     q_warning("\nWARNING(#%d):The max number of bad iterations exceeded. Can't decompose the current time series\n", skipCurrentPixel);
			     break;
			}

		}
		/*********************************/
		// WHILE(chainNumber<chainNumber)
		/*********************************/

		__END_IF_NOT_SKIP_TIMESESIRIES__  
	}

	/******************************************************/
	//
	// Finish all the chains and now acverage all of them
	//
	/******************************************************/

		// Average the results from multiple chains
		if (MCMC_CHAINNUM >= 2 && !skipCurrentPixel)
		{
			I32  N = opt->io.N;			
	
			I32   sMAXNUMKNOT = MODEL.sid >= 0 || MODEL.vid >= 0 ? MODEL.b[0].prior.maxKnotNum : -9999999;
			I32   tMAXNUMKNOT = MODEL.tid >= 0 ? MODEL.b[MODEL.tid].prior.maxKnotNum : -9999999;
			I32   oMAXNUMKNOT = MODEL.oid >= 0 ? MODEL.b[MODEL.oid].prior.maxKnotNum : -9999999;

			const F32 invChainNumber = 1.f /(F32)MCMC_CHAINNUM;

			#define _1(x)      *((F32PTR)result.x)*=invChainNumber
			#define _N(x)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N)
			#define _Nq(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,N*q)
			#define _q(x)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,q)
			#define _q2(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,q*q)
			#define _2N(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N+N)
			#define _2Nq(x)    r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N*q+N*q)
			#define _skn_1(x)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, sMAXNUMKNOT + 1)
			#define _tkn_1(x)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, tMAXNUMKNOT + 1)
			#define _okn_1(x)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, oMAXNUMKNOT + 1)

			F32 maxncpProb;	 
			_1(marg_lik);			
			_q2(sig2);
			if (MODEL.sid >= 0 || MODEL.vid>0) {
				_1(sncp); _skn_1(sncpPr);	     _N(scpOccPr); _Nq(sY); _Nq(sSD); 
				if (extra.computeSeasonOrder)   { _N(sorder); }
				if (extra.computeSeasonAmp)     {_N(samp), _N(sampSD);}
				if (extra.computeCredible)      { _2Nq(sCI); }
 
				*result.sncp_mode   = f32_maxidx(result.sncpPr,       sMAXNUMKNOT + 1, &maxncpProb);
				*result.sncp_median = GetPercentileNcp(result.sncpPr, sMAXNUMKNOT + 1, 0.5);
				*result.sncp_pct90  = GetPercentileNcp(result.sncpPr, sMAXNUMKNOT + 1, 0.9);
				*result.sncp_pct10  = GetPercentileNcp(result.sncpPr, sMAXNUMKNOT + 1, 0.1);
			}

			if (MODEL.tid >= 0) {
				_1(tncp); _tkn_1(tncpPr);	     _N(tcpOccPr); _Nq(tY); _Nq(tSD); 
				if (extra.computeTrendOrder)    { _N(torder); }
				if (extra.computeTrendSlope)    { _N(tslp), _N(tslpSD), _N(tslpSgnPosPr), _N(tslpSgnZeroPr);}
				if (extra.computeCredible)      { _2Nq(tCI); }

				*result.tncp_mode   = f32_maxidx(result.tncpPr,       tMAXNUMKNOT + 1, &maxncpProb);
				*result.tncp_median = GetPercentileNcp(result.tncpPr, tMAXNUMKNOT + 1, 0.5);
				*result.tncp_pct90  = GetPercentileNcp(result.tncpPr, tMAXNUMKNOT + 1, 0.9);
				*result.tncp_pct10  = GetPercentileNcp(result.tncpPr, tMAXNUMKNOT + 1, 0.1);
			}

			if (MODEL.oid >= 0) {
				_1(oncp); _okn_1(oncpPr);	   
				_N(ocpOccPr); _Nq(oY); _Nq(oSD);
				if (extra.computeCredible)      _2Nq(oCI);

				*result.oncp_mode   = f32_maxidx(result.oncpPr,       oMAXNUMKNOT + 1, &maxncpProb);
				*result.oncp_median = GetPercentileNcp(result.oncpPr, oMAXNUMKNOT + 1, 0.5);
				*result.oncp_pct90  = GetPercentileNcp(result.oncpPr, oMAXNUMKNOT + 1, 0.9);
				*result.oncp_pct10  = GetPercentileNcp(result.oncpPr, oMAXNUMKNOT + 1, 0.1);
			}

			/****************************************************************/

			if (extra.tallyPosNegSeasonJump && MODEL.sid >= 0) {
				_1(spos_ncp);             _1(sneg_ncp); //NEWLY ADDED					
				_skn_1(spos_ncpPr);    _skn_1(sneg_ncpPr); //NEWLY ADDED
				_N(spos_cpOccPr);         _N(sneg_cpOccPr); //NEWLY ADDED	
			}


			if (extra.tallyPosNegTrendJump) {
				_1(tpos_ncp);            _1(tneg_ncp); //NEWLY ADDED					
				_tkn_1(tpos_ncpPr); _tkn_1(tneg_ncpPr); //NEWLY ADDED
				_N(tpos_cpOccPr);         _N(tneg_cpOccPr); //NEWLY ADDED 
			}


			if (extra.tallyIncDecTrendJump) {
				_1(tinc_ncp);         _1(tdec_ncp); //NEWLY ADDED					
				_tkn_1(tinc_ncpPr); _tkn_1(tdec_ncpPr); //NEWLY ADDED
				_N(tinc_cpOccPr);     _N(tdec_cpOccPr); //NEWLY ADDED
			}

			if (extra.tallyPosNegOutliers && MODEL.oid >= 0) {
				_1(opos_ncp);            _1(oneg_ncp); //NEWLY ADDED					
				_okn_1(opos_ncpPr); _okn_1(oneg_ncpPr); //NEWLY ADDED
				_N(opos_cpOccPr);         _N(oneg_cpOccPr); //NEWLY ADDED 
			}

 
			#undef _1
			#undef _N
			#undef _2N
			#undef _skn_1
			#undef _tkn_1
			#undef _okn_1
		}
		if (!skipCurrentPixel) {		
			I32  N = opt->io.N;  
			if (MODEL.sid >= 0||MODEL.vid>0) 				tsRemoveNaNs(result.sSD, N);
			if (MODEL.tid >= 0)								tsRemoveNaNs(result.tSD, N);
			if (MODEL.tid >= 0 && extra.computeTrendSlope)  tsRemoveNaNs(result.tslpSD, N);		 
			if (MODEL.oid >= 0) 							tsRemoveNaNs(result.oSD, N);
		}

		// Compute Changepoints and their confidence intervals
		if (!skipCurrentPixel)
		{
			I32     N         = opt->io.N;
			F32     nan       = getNaN();     //A sloppy way to get a NAN

			F32  	threshold = 0.01f;
			F32PTR	mem       = Xnewterm;  //Xnewterm must have a length larger than or equla to 5*N + max(sncp, tncp);
			I32PTR  cptList   = (I32PTR)mem + 5LL * N;
			F32PTR  cptCIList = (F32PTR)mem + 6LL * N;

			I32   cptNumber;
			I32   trueCptNumber;
			F32   maxncpProb;
			const F32 T0 = (F32)opt->io.meta.startTime;
			const F32 dT = (F32)opt->io.meta.deltaTime;

			I32   sMAXNUMKNOT = MODEL.sid >= 0 || MODEL.vid >= 0 ? MODEL.b[0].prior.maxKnotNum : -9999999;
			I32   tMAXNUMKNOT = MODEL.tid >= 0 ? MODEL.b[MODEL.tid].prior.maxKnotNum : -9999999;
			I32   oMAXNUMKNOT = MODEL.oid >= 0 ? MODEL.b[MODEL.oid].prior.maxKnotNum : -9999999;

			I32   sMINSEPDIST = MODEL.sid >= 0 || MODEL.vid >= 0 ? MODEL.b[0].prior.minSepDist : -9999999;
			I32   tMINSEPDIST = MODEL.tid >= 0                   ? MODEL.b[MODEL.tid].prior.minSepDist : -9999999;
			I32   oMINSEPDIST = MODEL.oid >= 0                   ? 0 : -9999999;  //changed from 1 to 0 to pick up consecutive outlier spikes.

			I32   sLeftMargin = MODEL.sid >= 0 || MODEL.vid >= 0 ? MODEL.b[0].prior.leftMargin : -9999999;
			I32   tLeftMargin = MODEL.tid >= 0                   ? MODEL.b[MODEL.tid].prior.leftMargin : -9999999;
			I32   oLeftMargin = MODEL.oid >= 0                   ? 0 : -9999999;  //changed from 1 to 0 to pick up consecutive outlier spikes.


			I32   sRightMargin = MODEL.sid >= 0 || MODEL.vid >= 0 ? MODEL.b[0].prior.rightMargin : -9999999;
			I32   tRightMargin = MODEL.tid >= 0                   ? MODEL.b[MODEL.tid].prior.rightMargin : -9999999;
			I32   oRightMargin = MODEL.oid >= 0                   ? 0 : -9999999;  //changed from 1 to 0 to pick up consecutive outlier spikes.

			//--------------Season----------------------------
			if (extra.computeSeasonChngpt && (MODEL.sid >=0 || MODEL.vid>0 ))  	{
				cptNumber     = sMAXNUMKNOT;
				trueCptNumber = FindChangepoint_LeftRightMargins((F32PTR)result.scpOccPr, mem, threshold, cptList, cptCIList, N, sMINSEPDIST, cptNumber,sLeftMargin, sRightMargin);
				// In the returned result, cptList is a list of detected changepoints, which are all zero-based indices.
				// That is, a changepoint occurring at t1 has a value of 0.

				for (int i = 0; i < trueCptNumber; i++) {
					*(result.scp + i)	  = (F32)(*(cptList + i)) * dT + T0;
					*(result.scpPr + i)   = (F32)mem[i];
					I32 cptLoc = cptList[i] == 0 ? 1 : cptList[i];
					for (int j = 0; j < q; ++j) {
						*(result.scpAbruptChange +j*sMAXNUMKNOT + i) =	result.sY[j*N+cptLoc] - result.sY[j*N+cptLoc - 1];
					}					
				}
				for (int i = trueCptNumber; i < sMAXNUMKNOT; i++) {
					*(result.scp + i)             = nan,
					*(result.scpPr + i)           = nan;
					for (int j = 0; j < q; ++j) {
						*(result.scpAbruptChange + j * sMAXNUMKNOT+ i) = nan;
					}					
				}
				for (int i = 0; i < trueCptNumber; i++)
					*(result.scpCI + i) = (F32)(*(cptCIList + i)) * dT + T0,
					*(result.scpCI + sMAXNUMKNOT + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;
				for (int i = trueCptNumber; i < sMAXNUMKNOT; i++)
					*(result.scpCI + i) = nan,
					*(result.scpCI + sMAXNUMKNOT + i) = nan;
			}

			//--------------Trend----------------------------
			if (extra.computeTrendChngpt) {
				cptNumber     = tMAXNUMKNOT;
				trueCptNumber = FindChangepoint_LeftRightMargins((F32PTR)result.tcpOccPr, mem, threshold, cptList, cptCIList, N, tMINSEPDIST, cptNumber, tLeftMargin, tRightMargin);
				for (int i = 0; i < trueCptNumber; i++) {
					*(result.tcp + i)          = (F32)(*(cptList + i)) * dT + T0,
					*(result.tcpPr + i)         = (F32)mem[i];
					I32 cptLoc = cptList[i] == 0 ? 1 : cptList[i];
					for (int j = 0; j < q; ++j) {
						*(result.tcpAbruptChange + j*tMAXNUMKNOT + i) = result.tY[j * N + cptLoc] - result.tY[j * N + cptLoc - 1];
					}
				}
				for (int i = trueCptNumber; i < tMAXNUMKNOT; i++) {
					*(result.tcp + i)   = nan,
					*(result.tcpPr + i) = nan;
					for (int j = 0; j < q; ++j) {
						*(result.tcpAbruptChange + j * tMAXNUMKNOT + i) = nan;
					}					
				}					
				for (int i = 0; i < trueCptNumber; i++)
					*(result.tcpCI + i) = (F32)(*(cptCIList + i)) * dT + T0,
					*(result.tcpCI + tMAXNUMKNOT + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;
				for (int i = trueCptNumber; i < tMAXNUMKNOT; i++)
					*(result.tcpCI + i) = nan,
					*(result.tcpCI + tMAXNUMKNOT + i) = nan;
			}
	 

			/**************************************************************************************************/
			#define GET_CHANGPOINTS(NcpProb, KNOTNUM, MINSEP, LeftMargin, RightMargin, MAX_KNOTNUM, Y, CpOccPr, CP, CPPROB, CP_CHANGE, CP_CI)    \
			cptNumber     = MAX_KNOTNUM;  \
			trueCptNumber = FindChangepoint_LeftRightMargins((F32PTR)CpOccPr, mem, threshold, cptList, cptCIList, N, MINSEP, cptNumber, LeftMargin, RightMargin);\
			for (int i = 0; i < trueCptNumber; i++) {\
				*(CP + i)        = (F32) cptList[i]* dT + T0,\
				*(CPPROB+ i)     = (F32) mem[i];\
		         I32 cptLoc      = cptList[i] == 0 ? 1 : cptList[i];\
				 *(CP_CHANGE + i) = Y[cptLoc] - Y[cptLoc - 1];\
			}\
			for (int i = trueCptNumber; i <MAX_KNOTNUM; i++) {\
				*(CP        + i) = nan;\
				*(CPPROB    + i) = nan;\
				  *(CP_CHANGE + i) = nan;\
			}	\
			for (int i = 0; i < trueCptNumber; i++)\
				*(CP_CI+ i)                = (F32)cptCIList[i] * dT + T0,\
				*(CP_CI + MAX_KNOTNUM + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;\
			for (int i = trueCptNumber; i < MAX_KNOTNUM; i++)\
				*(CP_CI + i)               = nan,\
				*(CP_CI + MAX_KNOTNUM + i) = nan;
			/**************************************************************************************************/

			if (extra.tallyPosNegSeasonJump && MODEL.sid >= 0) {
				GET_CHANGPOINTS(result.spos_ncpPr ,result.spos_ncp, sMINSEPDIST, sLeftMargin, sRightMargin, sMAXNUMKNOT, result.samp,
					result.spos_cpOccPr, result.spos_cp, result.spos_cpPr, result.spos_cpAbruptChange, result.spos_cpCI);

				GET_CHANGPOINTS(result.sneg_ncpPr,result.sneg_ncp, sMINSEPDIST, sLeftMargin, sRightMargin, sMAXNUMKNOT, result.samp,
					result.sneg_cpOccPr, result.sneg_cp, result.sneg_cpPr, result.sneg_cpAbruptChange, result.sneg_cpCI);
			}			
			if (extra.tallyPosNegTrendJump) {
				GET_CHANGPOINTS(result.tpos_ncpPr, result.tpos_ncp, tMINSEPDIST, tLeftMargin, tRightMargin, tMAXNUMKNOT, result.tY,
					result.tpos_cpOccPr, result.tpos_cp, result.tpos_cpPr, result.tpos_cpAbruptChange, result.tpos_cpCI);

				GET_CHANGPOINTS(result.tneg_ncpPr, result.tneg_ncp, tMINSEPDIST, tLeftMargin, tRightMargin, tMAXNUMKNOT, result.tY,
					result.tneg_cpOccPr, result.tneg_cp, result.tneg_cpPr, result.tneg_cpAbruptChange, result.tneg_cpCI);
			}
			if (extra.tallyIncDecTrendJump) {
				GET_CHANGPOINTS(result.tinc_ncpPr,result.tinc_ncp, tMINSEPDIST, tLeftMargin, tRightMargin, tMAXNUMKNOT, result.tslp,
					result.tinc_cpOccPr, result.tinc_cp, result.tinc_cpPr, result.tinc_cpAbruptChange, result.tinc_cpCI);

				GET_CHANGPOINTS(result.tdec_ncpPr,result.tdec_ncp, tMINSEPDIST, tLeftMargin, tRightMargin, tMAXNUMKNOT, result.tslp,
					result.tdec_cpOccPr, result.tdec_cp, result.tdec_cpPr, result.tdec_cpAbruptChange, result.tdec_cpCI);
			}

		
			/**************************************************************************************************/
			#define OGET_CHANGPOINTS(NcpProb, KNOTNUM,MINSEP,LeftMargin, RightMargin, MAX_KNOTNUM, PROBCURVE, CP, CPPROB,CP_CI)    \
			cptNumber     = MAX_KNOTNUM;  \
			trueCptNumber = FindChangepoint_LeftRightMargins((F32PTR)PROBCURVE, mem, threshold, cptList, cptCIList, N, MINSEP, cptNumber, LeftMargin, RightMargin);\
			for (int i = 0; i < trueCptNumber; i++) {\
				*(CP + i)        = (F32) cptList[i]* dT + T0;\
				*(CPPROB+ i)     = (F32) mem[i];\
		         I32 cptLoc      = cptList[i] == 0 ? 1 : cptList[i];\
			}\
			for (int i = trueCptNumber; i <MAX_KNOTNUM; i++) {\
				*(CP        + i) = nan;\
				*(CPPROB    + i) = nan;\
		    }\
			for (int i = 0; i < trueCptNumber; i++) \
				*(CP_CI+ i)                = (F32)cptCIList[i] * dT + T0,\
				*(CP_CI + MAX_KNOTNUM + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;\
			for (int i = trueCptNumber; i < MAX_KNOTNUM; i++)\
				*(CP_CI + i)               = nan,\
				*(CP_CI + MAX_KNOTNUM + i) = nan;
			/**************************************************************************************************/

			if (extra.computeOutlierChngpt) {
				OGET_CHANGPOINTS(result.oncpPr,result.oncp, oMINSEPDIST, oLeftMargin, oRightMargin, oMAXNUMKNOT, result.ocpOccPr, result.ocp, result.ocpPr, result.ocpCI);
			}

			if (extra.tallyPosNegOutliers && MODEL.oid >=0) {
				OGET_CHANGPOINTS(result.opos_ncpPr,result.opos_ncp, oMINSEPDIST, oLeftMargin, oRightMargin, oMAXNUMKNOT,
					result.opos_cpOccPr, result.opos_cp, result.opos_cpPr, result.opos_cpCI);

				OGET_CHANGPOINTS(result.oneg_ncpPr,result.oneg_ncp, oMINSEPDIST, oLeftMargin, oRightMargin, oMAXNUMKNOT,
					result.oneg_cpOccPr, result.oneg_cp, result.oneg_cpPr, result.oneg_cpCI);
			}


		}
		

		// Recover the orignal Y values if the pixel is not skipped; the value will be needed below 
		// to dump into the output and compute R2 and RMSE. Here Xnewterm is used a buff
		// should not be touched until after the computation of R2 and RMSE
		if ( !skipCurrentPixel) {	
			I32  N  = opt->io.N;	
			I32  Nq = N * q;  // For MRBEAST

			for (int i = 0; i < q; ++i) {
				memcpy(Xnewterm+N*i, yInfo.Y + N*i, sizeof(F32)* N);
				f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], Xnewterm+N*i, N);

				F32PTR NULL_BUF_FOR_VALUES = (Xnewterm + N * i) + N;
				if (yInfo.nMissing > 0) {
					f32_gatherVec_scatterVal_byindex(Xnewterm + N*i, yInfo.rowsMissing, NULL_BUF_FOR_VALUES, getNaN(), yInfo.nMissing);
				}
			}			
		}

		if (opt->extra.dumpInputData && result.data != NULL) {
			I32  N = opt->io.N;
			I32  Nq = N * q;  // For MRBEAST
			memcpy(result.data, Xnewterm, sizeof(F32)* Nq);
		}

		if (opt->extra.dumpMCMCSamples) {
			result.smcmc = opt->io.out.result->smcmc;
			result.tmcmc = opt->io.out.result->tmcmc;
			result.omcmc = opt->io.out.result->omcmc;

			I32  N = opt->io.N;
			for (I32 i = 0; i < q; i++) {
				if (result.smcmc) {
					f32_mul_val_inplace(yInfo.sd[i], result.smcmc, N* MCMC_SAMPLES* MCMC_CHAINNUM);
					result.smcmc += N * MCMC_SAMPLES * MCMC_CHAINNUM;
				}				
				if (result.omcmc) {
					f32_mul_val_inplace(yInfo.sd[i], result.omcmc, N* MCMC_SAMPLES* MCMC_CHAINNUM);
					result.omcmc += N * MCMC_SAMPLES * MCMC_CHAINNUM;
				}				

				f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], result.tmcmc, N * MCMC_SAMPLES * MCMC_CHAINNUM);
				result.tmcmc += N * MCMC_SAMPLES * MCMC_CHAINNUM;
			}		

			if (opt->io.out.dtype == DATA_DOUBLE) {	
				result.smcmc = opt->io.out.result->smcmc;
				result.tmcmc = opt->io.out.result->tmcmc;
				result.omcmc = opt->io.out.result->omcmc;
				if (result.smcmc) f32_to_f64_inplace(result.smcmc, N * MCMC_SAMPLES * MCMC_CHAINNUM * q);
				if (result.omcmc) f32_to_f64_inplace(result.smcmc, N * MCMC_SAMPLES * MCMC_CHAINNUM * q);
				f32_to_f64_inplace(result.tmcmc, N* MCMC_SAMPLES* MCMC_CHAINNUM* q);
			}


		}

		// Compute R2 and RMSE
		// At this point, yInfo.Y is not used any longer, so is re-used
		// to stote the sume of fitted S,T and/or O.
		// Xnewterm still contains the orignal data and shiuldn't not be touched
  		if (!skipCurrentPixel) { 
			
			I32  N  = opt->io.N ;
			I32  Nq = N * q;  // For MRBEAST

			I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
			I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
			I08 hasTrendCmpnt   = 1;
			//Xnewterm is used at this point and should not be touched!
			F32PTR BUF      = yInfo.Y;
			f32_fill_val(0., BUF, Nq);
	
			if (hasTrendCmpnt)   f32_add_vec_inplace(result.tY, BUF, Nq);
			if (hasSeasonCmpnt)  f32_add_vec_inplace(result.sY, BUF, Nq);
			if (hasOutlierCmpnt) f32_add_vec_inplace(result.oY, BUF, Nq);
		
			for (int j = 0; j < q; ++j) {
				//For MRBEAST
				F32 r = f32_corr_rmse_nan(BUF+N*j, Xnewterm + N*j, N, &result.RMSE[j]);
				result.R2[j] = r * r;
			}
			
		}

		// Smooth the changepoint occurance probability curve
		if (!skipCurrentPixel) {
			I32  N = opt->io.N;

			I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
			I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
			I08 hasTrendCmpnt   = 1;
 
			F32PTR BUF = Xnewterm;
			if (hasTrendCmpnt && opt->extra.smoothCpOccPrCurve) {	
				memcpy(BUF, result.tcpOccPr, sizeof(F32)* N);
				f32_sumfilter(BUF, result.tcpOccPr,N, opt->prior.trendMinSepDist);
			}
			if (hasSeasonCmpnt && opt->extra.smoothCpOccPrCurve) {
				memcpy(BUF, result.scpOccPr, sizeof(F32)* N);
				f32_sumfilter(BUF, result.scpOccPr, N, opt->prior.seasonMinSepDist);
			}	
		}


	    if (skipCurrentPixel) BEAST2_Result_FillMEM(&result, opt, getNaN());
		/*********************************************/
		//Write outputs to the mem array or files	
		/*********************************************/	
		if (!skipCurrentPixel) {
			int N = opt->io.N;
			for (int i = 0; i < q; ++i) {
				if (yInfo.Yseason) {
				//If Y has been deseaonalized, add it back
					r_ippsAdd_32f_I(yInfo.Yseason + N*i,   result.sY+N* i, N);
					if (result.sCI) {
						r_ippsAdd_32f_I(yInfo.Yseason + N * i, result.sCI + 2*N*i, N);
						r_ippsAdd_32f_I(yInfo.Yseason + N * i, result.sCI + 2*N*i+ N, N);
					}
					if (extra.dumpInputData)
						r_ippsAdd_32f_I(yInfo.Yseason + N * i, result.data+N*i, N);
				}
				if (yInfo.Ytrend) {
				//If Y has been detrended, add it back
					r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.tY + N*i, N);
					if (result.tCI) {
						r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.tCI+ 2*N*i,   N);
						r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.tCI+ 2*N*i+N, N);
					}
					if (extra.dumpInputData)
						r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.data + N * i, N);
				}
			} //for (int i = 0; i < q; ++i)
		}

	 
		BEAST2_WriteOutput(opt, &result, pixelIndex);
		 
		  

		//if (!skipCurrentPixel)	NUM_OF_PROCESSED_GOOD_PIXELS++; //avoid the branch
		NUM_OF_PROCESSED_GOOD_PIXELS += !skipCurrentPixel;              //this is a global variable.
		NUM_OF_PROCESSED_PIXELS++;					//this is also a global variable.


		F64 elaspedTime = GetElaspedTimeFromBreakPoint();
		if (NUM_OF_PROCESSED_GOOD_PIXELS > 0 && NUM_PIXELS > 1 && (pixelIndex % 50 == 0 || elaspedTime > 1)) 		{
			F64 estTimeForCompletion = GetElapsedSecondsSinceStart()/NUM_OF_PROCESSED_GOOD_PIXELS * (NUM_PIXELS - pixelIndex);
			printProgress2((F32)pixelIndex / NUM_PIXELS, estTimeForCompletion, extra.consoleWidth, Xnewterm, 0);
			if (elaspedTime > 1) SetBreakPointForStartedTimer();
		}

		#if DEBUG_MODE == 1
		r_printf("TREND: birth%4d/%-5d|death%4d/%-5d|merge%4d/%-5d|move%4d/%-5d|chorder%4d/%-5d\n", 
			      accT[0], flagT[0] , accT[1], flagT[1], accT[2], flagT[2], accT[3], flagT[3], accT[4], flagT[4]);
		r_printf("SEASN: birth%4d/%-5d|death%4d/%-5d|merge%4d/%-5d|move%4d/%-5d|chorder%4d/%-5d\n",
			      accS[0], flagS[0], accS[1], flagS[1], accS[2], flagS[2], accS[3], flagS[3], accS[4], flagS[4]);
		#endif

	} //for (U32 pixelIndex = 1; pixelIndex <= TOTALNUMPIXELS; pixelIndex++)





	/***********************************************************/
	// This is the ending bracekt of the iteration through pixels
	/***********************************************************/

	r_vslDeleteStream(&stream);
	MEM.free_all(&MEM);
 
	return 1;
} /* End of beastST() */


#include "abc_000_warning.h"