#include "mex.h"

#if !defined(MAX)
#define MAX(A, B)       ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#endif

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]
		 )
{

    double *ExcG, *InhG, *Parameters, *Spikes, *Voltage, *CurrentNoise, *ptr;
    unsigned int Length, SpikeCount, t;
    double VRevInh, VRevExc, VRevLeak, GLeak, Cap, TStep, Threshold, AbsRefractTime, AHPDecay, AHPAmp, AHPVRev;
    double AHP;
    
    // check for proper number of arguments   
	if (nrhs != 5 || nlhs != 2){
	    mexErrMsgTxt("Usage: [Spikes,Voltage] = PredictRGCVoltageAndSpikes_c(ExcG, InhG, CurrentNoise, Parameters, Length\n");
		return;
	}
	
	// get number of data points for buffers
	ExcG = mxGetPr(prhs[0]);
	InhG = mxGetPr(prhs[1]);
    CurrentNoise = mxGetPr(prhs[2]);
	Parameters = mxGetPr(prhs[3]);
    Length = mxGetScalar(prhs[4]);
    
    if ((Spikes = (double *) malloc(sizeof(double) * Length)) == -1){
  	    mexErrMsgTxt("Malloc error, Spikes\n");
		return;
    }
    
    if ((Voltage = (double *) malloc(sizeof(double) * Length)) == -1){
  	    mexErrMsgTxt("Malloc error, Spikes\n");
		return;
    }
    
    VRevInh = Parameters[0];
    VRevExc = Parameters[1];
    VRevLeak = Parameters[2];
    GLeak = Parameters[3];
    Cap = Parameters[4];
    TStep = Parameters[5];
    Threshold = Parameters[6];
    AbsRefractTime = Parameters[7];
    AHPDecay = Parameters[8];
    AHPAmp = Parameters[9];
    AHPVRev = Parameters[10];

    Voltage[0] = VRevLeak;
    
    AHP = 0;
    SpikeCount = 0;
    for (t = 1; t < Length; t++){
        Voltage[t] = Voltage[t-1] + TStep * (CurrentNoise[t-1] - AHP * (Voltage[t-1] - AHPVRev) - ExcG[t-1] * (Voltage[t-1] - VRevExc) - InhG[t-1] * (Voltage[t-1] - VRevInh) - GLeak * (Voltage[t-1] - VRevLeak)) / Cap;
        AHP = AHP - TStep * AHP / AHPDecay;
        if (Voltage[t] > Threshold){
            if (SpikeCount == 0){
                Spikes[0] = t * TStep;
                SpikeCount += 1;    
                AHP = AHP + AHPAmp;
            } else{ 
                if ((t*TStep - Spikes[SpikeCount-1]) > AbsRefractTime){
                    Spikes[SpikeCount] = t * TStep;
                    SpikeCount += 1;
                    AHP += AHPAmp;
                }
            }
        }
    }

    plhs[0] = mxCreateDoubleMatrix(1, SpikeCount, mxREAL);
	ptr = mxGetPr(plhs[0]);

    memcpy(ptr, Spikes, sizeof(double)*SpikeCount);
 
    plhs[1] = mxCreateDoubleMatrix(1, Length, mxREAL);
	ptr = mxGetPr(plhs[1]);

    memcpy(ptr, Voltage, sizeof(double)*Length);

    return;

}
