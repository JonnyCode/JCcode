/* 
   Matlab call: [V F E U M] = gpn_ada_mex (D, B, G, ep, dt, U_init, m_init)
   Description: Goal Programming Network with adaptive learning rate +
		exclusive constraints
		G = [nc x 1]
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "mex.h"
#include "matrix.h"

/* Global Variables */
double _tmp;

/* Macro Definitions */
#define sqr(x) ((_tmp=(x))*_tmp)

/* Function Prototypes */
double std (double*, int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nc, ns, ind1, ind2, e, ep, k, n;
	double dt, m, m_max, m_min, m_init, st_conv;
	double *ARG, *D, *B, *G, *U_init, *U, *V, *E, *M,
	       *F, *GF, *V_out, *conv_flag;
	double *AUX1, *AUX2;

	/**********************************/
	/*   Check number of arguments.   */
	/**********************************/
	if(nrhs != 7)
		mexErrMsgTxt("Wrong number of input arguments");
	else if(nlhs>6)
		mexErrMsgTxt("Too many output arguments!");

	/**********************************/
	/*       Initialise Inputs        */
	/**********************************/
	ind1 = -1;
	D = mxGetPr(prhs[++ind1]);
	nc = (int) mxGetM(prhs[ind1]);
	ns = (int) mxGetN(prhs[ind1]);

	B = mxGetPr(prhs[++ind1]);
	G = mxGetPr(prhs[++ind1]);
	ARG = mxGetPr(prhs[++ind1]);
	ep = (int) *ARG;
	ARG = mxGetPr(prhs[++ind1]);
	dt = *ARG;
	U_init = mxGetPr(prhs[++ind1]);
	ARG = mxGetPr(prhs[++ind1]);
	m_init = *ARG;

	
	/* ****************************** */
	/*       Initialise Outputs       */
	/* ****************************** */
	ind1 = -1;
	plhs[++ind1] = mxCreateDoubleMatrix(1,ns, mxREAL);
	V_out = (double*) mxGetPr(plhs[ind1]);
	plhs[++ind1] = mxCreateDoubleMatrix(1,nc, mxREAL);
	F = (double*) mxGetPr(plhs[ind1]);
	plhs[++ind1] = mxCreateDoubleMatrix(ep,1, mxREAL);
	E = (double*) mxGetPr(plhs[ind1]);
	plhs[++ind1] = mxCreateDoubleMatrix(ns,ep+1, mxREAL);
	U = (double*) mxGetPr(plhs[ind1]);
	plhs[++ind1] = mxCreateDoubleMatrix(ep,1, mxREAL);
	M = (double*) mxGetPr(plhs[ind1]);
	plhs[++ind1] = mxCreateDoubleMatrix(1,1, mxREAL);
	conv_flag = (double*) mxGetPr(plhs[ind1]);
	conv_flag[0] = .0;

	/* ******************** */
	/* Array Initialisation */
	/* ******************** */
	GF = (double*) mxCalloc (nc, sizeof(double));

	AUX1 = (double*) mxCalloc (ns, sizeof(double));
	AUX2 = (double*) mxCalloc (nc, sizeof(double));

	/* ******************** */
	/* Determine Parameters */
	/* ******************** */
	m_max = 25;
	m_min = .1;

	/**********************************/
	/*            Algorithm           */
	/**********************************/

			/* Initial Training */
	m = m_init;
	M[0] = m_init;
	e = 1;
	for (k=0; k<ns; k++)
		U[k] = U_init[k];

	/* Compute Network States */
	V = &U[0];

	/* Output F amps */
	for (n=0; n<nc; n++)
		F[n] =- B[n];
	for (k=ind1=0; k<ns; k++)
		for (n=0; n<nc; n++,ind1++)
			F[n] += D[ind1]*V[k];
	for (n=E[e-1]=0; n<nc; n++)
		E[e-1] += fabs(F[n]);

	/* Update Network States */
	for (k=ind1=0; k<ns; k++)
		for (n=AUX1[k]=0; n<nc; n++,ind1++)
			AUX1[k] -= F[n]*D[ind1];

	for (k=0; k<ns; k++)
		U[e*ns+k] = U[(e-1)*ns+k] + (m*dt)*AUX1[k];
	M[e] = m;

			/* Actual Training */
	for (e=2,ind2=2; e<ep; e++,ind2+=ns)
		{
		/* printf ("*\n"); */
		/* Compute Network States */
		V = &U[ind2-ns];

		/* Output F amps */
		for (n=0; n<nc; n++)
			F[n] =- B[n];
		for (k=ind1=0; k<ns; k++)
			for (n=0; n<nc; n++,ind1++)
				F[n] += D[ind1]*V[k];
		for (n=E[e-1]=0; n<nc; n++)
			{
			GF[n] = F[n]*G[n];
			E[e-1] += fabs(GF[n]);
			}

		/* Adaptive Learning Rate */
		if ((E[e-2]/E[e-1])>.9999)
			{
			m *= 1.05;
			if (m>m_max)
				m = m_max;

			/* Update Network States */
			for (k=ind1=0; k<ns; k++)
				for (n=AUX1[k]=0; n<nc; n++,ind1++)
					AUX1[k] -= GF[n]*D[ind1];
			for (k=0; k<ns; k++)
				U[ind2+k] = U[ind2-ns+k] + (m*dt)*AUX1[k];
			}
		else
			{
			m /= 2.0;
			if (m<m_min)
				m = m_min;
			for (k=0; k<ns; k++)
				U[ind2+k] = U[ind2-ns+k];
			}
		M[e] = m;

		/* Check Convergence */
		if ((e>15) && (std(&E[e-11],10)<1e-10))
			{
			conv_flag[0] = 1.0;
			break;
			}
		}

	for (k=0; k<ns; k++)
		V_out[k] = V[k];

	mxFree(GF); mxFree (AUX1); mxFree(AUX2);
}


double std (double *X, int l)
{
	int k;
	double me=0,st=0;
	double *X2;

	X2 = X;
	/* for (k=0; k<l; me+=*(X2++),k++); */
	for (k=0; k<l; k++,X2++)
		me += *(X2);
	me /= (double) l;
	X2 = X;
	/* for (k=0; k<l; st+=sqr(*(X2++)-me),k++); */
	for (k=0; k<l; k++,X2++)
		st += sqr(*(X2)-me);
	st /= (double) (l-1);

	return (sqrt(st));

}

