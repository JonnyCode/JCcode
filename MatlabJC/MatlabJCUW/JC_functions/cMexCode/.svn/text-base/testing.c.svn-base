/*
 *  testing.c
 *  
 *
 *  Created by Jon Cafaro on 5/6/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double y,x,*q,z,u ;	/* declares y,x,z are double precision variables, and q is a pointer */
	mxArray *out, *structure;		/* declares out as an array */
	
	structure = mxGetField(prhs[2],0,"testStruct"); /* this gets the value from the field named testStruct and puts in mxArray called structure */
	
	y=mxGetScalar(prhs[0]);		/* y is a scalar taken from the parameter right-hand side (prhs- first input in matlab)  */
	x=mxGetScalar(prhs[1]);
	
	u = mxGetScalar(structure); /* get scalar from struture */
	
	out=mxCreateDoubleMatrix(1,1,mxREAL);	/* creates an array of memmory */
	q = mxGetPr(out);	/* q is a pointer to the memmory array out */
	
	z= x+y+u;	/* here we are making z the sum of two doubles */
	
	q[0]=z;	/* assigns z to the location of memmory specified by q (the q[0] is called "dereferencing" and is the same as *q) */
	
	plhs[0] = out; /* now we are making the output equal to the the stuff in the memmory "out" */
	
	if(4< z){				/* if loop */
		mexPrintf("x+y+u is bigger than 4\n");

	}
}
	


	
	
/*	mexPrintf("You gave me %i parameters.\n", nrhs);
	mexPrintf("First thing you gave me is an %s.", is);*/










