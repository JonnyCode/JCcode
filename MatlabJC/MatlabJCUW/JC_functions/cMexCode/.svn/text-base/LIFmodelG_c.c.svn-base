/*
 *  LIFmodelG_c.c
 *  
 *
 *  Created by Jon Cafaro on 5/8/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include <mex.h>
#include "matrix.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray  *current, *time, *voltage, *threshold ;
	double trial, *Gexc, *Ginh, *V, *I, ref, *Thresh ;
	int  NumTrials, NumSampPnts, t, x ;
	mxArray *Gleak_c, *Eleak_c, *Eexc_c, *Einh_c, *SampleRate_c, *Cap_c, *RefAmp_c, *RefDecay_c, *RefAbs_c, *Vthresh_c ;
	double Gleak, Eleak, Eexc, Einh, Eahp, SampleRate, Cap, RefAmp, RefDecay, RefAbs, Vthresh ;
	
	Gexc = mxGetPr(prhs[0]) ; /* get pointers to input matricies */
	Ginh = mxGetPr(prhs[1]) ;
	
	Gleak_c = mxGetField(prhs[2],0,"Gleak") ; /*get mx arrays from field strucutures */
	Eleak_c = mxGetField(prhs[2],0,"Eleak") ;
	Eexc_c = mxGetField(prhs[2],0,"Eexc") ;
	Einh_c = mxGetField(prhs[2],0,"Einh") ;
	SampleRate_c = mxGetField(prhs[2],0,"SampleRate") ;
	Cap_c = mxGetField(prhs[2],0,"cap") ;
	RefDecay_c = mxGetField(prhs[2],0,"RelRefTau") ;
	RefAmp_c = mxGetField(prhs[2],0,"RelRefAmp") ;
	RefAbs_c = mxGetField(prhs[2],0,"AbsRef") ;
	Vthresh_c = mxGetField(prhs[2],0,"Vthresh") ;

	Gleak = mxGetScalar(Gleak_c) ; /* get values from above mx arrays */
	Eleak = mxGetScalar(Eleak_c) ;
	Eexc = mxGetScalar(Eexc_c) ;
	Einh = mxGetScalar(Einh_c) ;
	SampleRate = mxGetScalar(SampleRate_c) ;
	Cap = mxGetScalar(Cap_c) ;
	RefDecay = mxGetScalar(RefDecay_c) ;
	RefAmp = mxGetScalar(RefAmp_c) ;
	RefAbs = mxGetScalar(RefAbs_c) ;
	Vthresh = mxGetScalar(Vthresh_c) ;
	
	NumTrials = mxGetM(prhs[0]) ; /* get number of trials from first matrix input */
	NumSampPnts = mxGetN(prhs[0]) ; /* get number of sample points from first matrix input*/
	
	if(NumTrials>1){
		 mexErrMsgTxt("this c code will only deal with one set of conductances at a time");
	}
		
	voltage = mxCreateDoubleMatrix(NumTrials,NumSampPnts,mxREAL) ;	/* creates an array of memmory */
	current = mxCreateDoubleMatrix(NumTrials,NumSampPnts,mxREAL) ;	/* creates an array of memmory */
	threshold = mxCreateDoubleMatrix(NumTrials,NumSampPnts,mxREAL) ;
	
	V = mxGetPr(voltage) ; /* get a pointer to that array of memmory */
	I = mxGetPr(current) ;
	Thresh = mxGetPr(threshold) ;
	
	for(t = 0 ; t<NumSampPnts ; t++){
		Thresh[t] = Vthresh ;
	}
	
	/* LIF code below */

	ref = 0 ; /* refractory period */
	V[0] = Eleak ; /* set rest voltage */
	I[0] = 0 ; /* there is no synaptic current (ignores first point in conductances) */
	
	for(t = 1 ; t<NumSampPnts ; t++){ /* for each time point */
		I[t] = (Gexc[t]*(V[t-1]-Eexc))+(Ginh[t]*(V[t-1]-Einh))+(Gleak*(V[t-1]-Eleak)) ; /* calculate synaptic current from previous voltage and current conductances */			
		if(ref==0){ /* if not withing the absolute refractory period */	
			V[t] = V[t-1] - I[t]/(Cap*SampleRate) ; /* calculate current voltage from previous voltage and current current */
		}
		else{
			V[t] = Eleak ;
			ref = ref - 1 ;
		}
		
		if(V[t]>Thresh[t]){ /* if above spike threshold */
			V[t] = 50 ; /* spike */
			ref = RefAbs * SampleRate ; /* absolute ref period */
			for(x = 0 ; x<(RefDecay*SampleRate*3) ; x++) /* for speed I have cut the decay off at longer times */
				if((t+x)>NumSampPnts)
				Thresh[t+x] = Thresh[t+x] + RefAmp * exp(x/(-1*RefDecay*SampleRate)) ; /* add relative refractory exponential decay to spike threshold */
			
		}
				
	}
		
	plhs[0] = voltage ; /* set voltage output matricies */
	plhs[1] = current ; /* set synaptic current output matricies */
	plhs[2] = threshold ;
}
