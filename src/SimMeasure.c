/****************************************************************************************************
* SimMeasure.c
*
* Authors: 	Shannon M. Bell
*          	Michigan State University
*          	Department of Biochemistry & Molecular Biology
*		Quantitative Biology
*
* Version 3.2	April 19, 2011
*
* Notes: SMB conceived of the algorithm, and wrote an early version in R. SMB used skelton code provided by
* Dr.Lyle Burgoon (EPA) to translate algorithm to C so that it can be run faster within R.
*******************************************************************************************************/



#include "Rdefines.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "math.h"

SEXP SimMeasure(SEXP data_matrix, SEXP thresh){
	// data_matrix is a matrix
	// thresh is a double

	int num_cols, num_rows;
	double *rx = REAL(data_matrix), *rans, t;
	SEXP retval;
	//d is defunct in this version for now
	//Check to make sure that everything is of the proper type before going further...
	if (isMatrix(data_matrix)){
		num_cols = ncols(data_matrix);
		num_rows = nrows(data_matrix);
	}
	else{
		REprintf("invalid matrix.\n");
		return R_NilValue;
	}
	if (isNull(thresh)){
		REprintf("warning, setting threshold to 0 by default.\n");
		t = 0;
	}
	else{
		t = REAL(thresh)[0];
	}

	//Check to see if the matrix has any null values
	if (isNull(data_matrix)){
		REprintf("matrix must not be NULL.\n");
		return R_NilValue;
	}

	PROTECT(retval = allocMatrix(REALSXP, num_cols, num_cols));
	rans = REAL(retval);

	for(int i = 0; i < num_cols; i++){
		for(int q = 0; q < num_cols; q++){
			double cor_val = 0.0;
			int count_row_nas = 0;
			for(int wi = 0; wi < num_rows; wi++){
				//REprintf("row_nas: %d\n", row_nas[wi]);
				//Check to see if we need to skip this row b/c BOTH of the elements are NA
				if(ISNAN(rx[i * num_rows + wi]) && ISNAN(rx[q * num_rows + wi])){
					count_row_nas++;
					continue;
				}
				//Check to see if we need to skip this row b/c BOTH of the elements are <hit
				//this is for cases where the 'non-hits' werent removed
				if(fabs(rx[i * num_rows + wi]) < t && fabs(rx[q * num_rows + wi]) < t){
					count_row_nas++;
					continue;
				}
			}
			//Calculate the parts of the similarity function
			double nm = 0.0; double pm = 0.0; double om = 0.0;
			double pcnt = 0.0; double ocnt = 0.0;
			for(int wi = 0; wi < num_rows; wi++){
				//Check to see if one of the elements are NA
				if(ISNAN(rx[i * num_rows + wi]) && ISNAN(rx[q * num_rows + wi])){
					continue;
				}
				else if(fabs(rx[i * num_rows + wi]) < t && fabs(rx[q * num_rows + wi]) < t){
					continue;
				}
				//one value missing, the other above threshold
				else if((ISNAN(rx[i * num_rows + wi]) && fabs(rx[q * num_rows + wi]) >= t) || (ISNAN(rx[q * num_rows + wi]) && fabs(rx[i * num_rows + wi]) >= t)){
					nm++;
				}
				//one value below threshold, other value above
				else if((fabs(rx[i * num_rows + wi]) < t && fabs(rx[q * num_rows + wi]) >= t) || (fabs(rx[q * num_rows + wi]) < t && fabs(rx[i * num_rows + wi]) >= t)){
					nm++;
				}
				//both are duds
				else if((ISNAN(rx[i * num_rows + wi]) && fabs(rx[q * num_rows + wi]) < t) || (ISNAN(rx[q * num_rows + wi]) && fabs(rx[i * num_rows + wi]) < t)){
					continue;
				}
				else{
					if(rx[i * num_rows + wi] * rx[q * num_rows + wi] >= 0){
						
						pm = pm + 1-(fabs(fabs(rx[i * num_rows + wi]) - fabs(rx[q * num_rows + wi]))/((fabs(rx[i * num_rows + wi]) + fabs(rx[q * num_rows + wi]))/2));
						pcnt++;
					}
					else{
						om = om + 1-(fabs(fabs(rx[i * num_rows + wi]) - fabs(rx[q * num_rows + wi]))/((fabs(rx[i * num_rows + wi]) + fabs(rx[q * num_rows + wi]))/2));
						ocnt++;
					}
				}
			}

			//Calculate the correlation, and the correlation matrix
			if(num_rows - count_row_nas < 1){
				cor_val = NA_REAL;
			}
			else{
				cor_val = (pm - om) / (pcnt + ocnt + nm/(pcnt+ocnt+nm));
			}
			rans[i * num_cols + q] = cor_val;
		}
	}
	UNPROTECT(1);
	return retval;				//return the correlation matrix
}

R_CallMethodDef callMethods[] =
{
    {"SimMeasure", (DL_FUNC)&SimMeasure, 1},
    {NULL,NULL, 0}
};

void R_init_pretty_matrix(DllInfo *dll)
{
    R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
}

