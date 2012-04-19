/****************************************************************************************************
* NetVal.c
*
* Authors: 	Shannon M. Bell
*          	Michigan State University
*          	Department of Biochemistry & Molecular Biology
*		Quantitative Biology
*
* Version 1	March 1, 2012
*
* Notes: SMB conceived of and wrote the algorithm.
* This algorithm is written to compare 2 adjacency matricies (correlation) and calculate
* the number of TP, FP, TN, FN based on communties defined from a heirctical tree
*******************************************************************************************************/



#include "Rdefines.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "math.h"

SEXP netVal(SEXP tree){
	// tree is a vector

	int l;
	double *rx = REAL(tree), *rans;
	l = length(tree);
	SEXP retval;
	PROTECT(retval = allocMatrix(REALSXP, l, l));
	rans = REAL(retval);
	for(int i = 0; i < l; i++){
		for(int j = 0; j < l; j++){
			if(rx[i] == rx[j]){
				rans[i * l + j] =1;
			}
			else{
				rans[i * l + j] =0;
			}
		}
	}
	UNPROTECT(1);
	return retval;				//return the correlation matrix
}
