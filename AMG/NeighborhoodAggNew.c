#include "mex.h"
#include <omp.h>
#include <math.h>


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex  -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims NeighborhoodAggNew.c
    mwIndex k,global_idx;
// // 
    mwIndex *C_S = mxGetIr(prhs[0]);
    mwIndex *starts_S = mxGetJc(prhs[0]);
    double  *V_S = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    mwIndex neighbor, chosen, agg_of_neighbor;
	double* aggr;
	double chosen_score = 0, avg_sparsity = 0;
	double* aux = (double*)malloc(n*sizeof(double));
	int* aux_count = (int*)malloc(n*sizeof(int));
    int Neighbors_Aggregated_flag = 0;
    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
	aggr = mxGetPr(plhs[0]);
	
    // initialization: aggr will hold the answer of the aggregation.
		
    for (k = 0 ; k < n ; k++){
        aggr[k] = 0;
		aux[k] = 0;
		aux_count[k] = 0;
		avg_sparsity += starts_S[k+1] - starts_S[k];
    }
	avg_sparsity /= n;
	// Here we set "too popular" variables to be ignored in the aggregation - so they won't be a seed.
	for (k = 0 ; k < n ; k++){
        if (starts_S[k+1] - starts_S[k] > 3*avg_sparsity){
			aux_count[k] = -1; 
		}
    }
    for (k = 0 ; k < n ; k++){
        Neighbors_Aggregated_flag = 0;
		if (aux_count[k] == -1){
			continue; // ignore...
		}
        for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
            if (aggr[C_S[global_idx]] != 0){
                Neighbors_Aggregated_flag = 1;
                break;
            }
        }
        if (Neighbors_Aggregated_flag==0){
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
				if (aux_count[C_S[global_idx]]!=-1){
					aggr[C_S[global_idx]] = k+1; //conversion to Matlab's indices;
					aux_count[k]++;
				}
            }
        }
    }
	// here we repeat the process except ignored variables...
	for (k = 0 ; k < n ; k++){
        Neighbors_Aggregated_flag = 0;
		if (aux_count[k]!=-1){
			continue;
		}
		aux_count[k] = 0;
        for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
            if (aggr[C_S[global_idx]] != 0){
                Neighbors_Aggregated_flag = 1;
                break;
            }
        }
        if (Neighbors_Aggregated_flag==0){
            for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
				aggr[C_S[global_idx]] = k+1; //conversion to Matlab's indices;
				aux_count[k]++;
            }
        }
    }
	
	for (k = 0 ; k < n ; k++){
		chosen_score = 0;
		if (aggr[k] == 0){
			for (global_idx = starts_S[k] ; global_idx < starts_S[k+1] ; global_idx++ ){
				if (aggr[C_S[global_idx]] > 0){ // aggr can be negative here....					
                    agg_of_neighbor = aggr[C_S[global_idx]]-1;
					aux[agg_of_neighbor] += V_S[global_idx];
                }
			}
			for (global_idx = starts_S[k] ; global_idx<starts_S[k+1] ; global_idx++ ){
				if (aggr[C_S[global_idx]]>0){ // aggr can be negative here....
                    agg_of_neighbor = aggr[C_S[global_idx]]-1;
                    if (chosen_score < aux[agg_of_neighbor]/aux_count[agg_of_neighbor]) {
						chosen_score = aux[agg_of_neighbor]/aux_count[agg_of_neighbor];
						chosen = agg_of_neighbor;
                        aux[agg_of_neighbor] = 0;
					}
                }
			}
            aggr[k] = -((double)(chosen+1)); // we don't want the new aggregate to spread - use only original aggregate
		}
	}
	for (k = 0 ; k < n ; k++){
		if (aggr[k] < 0){
			aggr[k] = -aggr[k];
		}
	}
    free(aux);
	free(aux_count);
}
 