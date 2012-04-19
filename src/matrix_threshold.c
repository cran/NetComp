/****************************************************************************************************
* matrix_threshold.c
* 
* Authors: 	Lyle D. Burgoon, Ph.D.  			and 		Shannon M. Bell
*
* Version 1.0	January 17, 2011	Initial Write
* Version 1.1	March 2, 2012		Expanded functionality
*
* This function will replace any values less than 'threshold' with the minval,
* those that are >= threshold are replaced with maxval. If maxval or minval = NULL
* then the origional value is left intact. if abs=TRUE then the absolute value
* of the observation is considered. If rm.na=FALSE then NA's are maintined, rm.na=TRUE will replace NA's with minval
* in the event that minval=NULL, na's will be replaced with 0.
* Copyright 2011-2012 Lyle D. Burgoon and Shannon M. Bell. All Rights Reserved.
******************************************************************************************************/

/******************************************************************************************************
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
********************************************************************************************************/

#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "Rdefines.h"
#include <stdlib.h>

/**
* matrix1 => adjacency matrix (must be square)
* matrix2 => adjacency matrix (must be square)

**/
SEXP matrix_threshold(SEXP data_matrix, SEXP threshold, SEXP minval, SEXP maxval, SEXP abs, SEXP rmna){
	int num_cols, num_rows, num_cells;
	double *rx = REAL(data_matrix), *rans, *threshold_value = REAL(threshold), mi, ma;
	Rboolean use_abs = asLogical(abs);
	Rboolean rm_nas = asLogical(rmna);
	if (isNull(minval) == FALSE){
		mi = REAL(minval)[0];
	}
	else{
		mi = 0; //has no effect on output
	}
	if (isNull(maxval) == FALSE){
		ma = REAL(maxval)[0];
	}
	else{
		ma = 1; //has no effect on output
	}
	//char *intarray = NULL;
	//check to see if matrix is valad and determine the total number of cells
	//for matrix search num_cells
	SEXP retval;
	if (isMatrix(data_matrix)){
		num_cols = ncols(data_matrix);
		num_rows = nrows(data_matrix);
		num_cells = num_cols * num_rows;
	}
	else{
		Rprintf("invalid matrix.\n");
		return R_NilValue;
	}
	
	//Check to see if the matrix has any null values
	if (isNull(data_matrix)){
		Rprintf("matrix must not be NULL.\n");
		return R_NilValue;
	}
	
	PROTECT(retval = allocMatrix(REALSXP, num_rows, num_cols));
	rans = REAL(retval);
	//printf("use_abs: %d\n", use_abs);
	//this loop goes through each cell in the input matrix
	//test to see if the value is below the threshold and replace acordingly
	for(int i = 0; i < num_cells; i++){
		//printf("col %d\n", i + 1);
		if(rm_nas){ //should NA values be removed and replaces with minval?
			if(use_abs){ //are we considering the absolute value of the cell?
				if((ISNAN(rx[i]) == TRUE || fabs(rx[i]) < threshold_value[0]) && isNull(minval) == FALSE){
					rx[i] = mi;
				}
				else if(fabs(rx[i]) >= threshold_value[0] && isNull(maxval) == FALSE){
					rx[i] = ma;
				}
				else if(ISNAN(rx[i]) == TRUE && isNull(minval) == TRUE){
					rx[i] = 0; //defaults to 0 if minval is not provided and rm.ma=TRUE
				}
				else {
					rx[i] = rx[i];
				}
			}
			else{
				if((ISNAN(rx[i]) == TRUE || rx[i] < threshold_value[0])&& isNull(minval) == FALSE){
					rx[i] = mi;
				}
				else if(rx[i] >= threshold_value[0] && isNull(maxval) == FALSE){
					rx[i] = ma;
				}
				else if(ISNAN(rx[i]) == TRUE && isNull(minval) == TRUE){
					rx[i] = 0;
				}
				else {
					rx[i] = rx[i];
				}
			}
		}
		else{ //if NA's are to be maintained
			if(use_abs){
				if(fabs(rx[i]) < threshold_value[0] && isNull(minval) == FALSE){
					rx[i] = mi;
				}
				else if(fabs(rx[i]) >= threshold_value[0] && isNull(maxval) == FALSE){
					rx[i] = ma;
				}
				else {
					rx[i] = rx[i];
				}
			}
			else{
				if(rx[i] < threshold_value[0]&& isNull(minval) == FALSE){
					rx[i] = mi;
				}
				else if(rx[i] >= threshold_value[0] && isNull(maxval) == FALSE){
					rx[i] = ma;
				}
				else {
					rx[i] = rx[i];
				}
			}
			
		}
		//printf("%g\n", threshold_value[0]);
		rans[i] = rx[i];
	}
	UNPROTECT(1);
	return retval;
}

//R_CallMethodDef callMethods[] =
//{
//    {"matrix_threshold", (DL_FUNC)&matrix_threshold, 1},
//    {NULL,NULL, 0}
//};
//
//void R_init_pretty_matrix(DllInfo *dll)
//{
//    R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
//}



