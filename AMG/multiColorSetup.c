#include "mex.h"

void swap(mwIndex* a, mwIndex* b){
    mwIndex t = *a;
    *a = *b;
    *b = t;
}
void makeLargestTwoColorsFirstAndLast(mwIndex* aux_colors, 
        mwIndex* aux_color_elems, mwIndex n_colors, mwIndex n );

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /*
     mex -largeArrayDims multiColorSetup.c
     */
    
    /*Input pointers and variables*/
    mwIndex i,j,k,n_colors = 0;
    mwIndex *aux_color_elems = 0;
    mwIndex *aux_colors = 0;
    
    mwIndex *C_t = mxGetIr(prhs[0]);
    mwIndex *starts_t = mxGetJc(prhs[0]);
    double* vals_t = mxGetPr(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    mwIndex maximalDegree = 0;
    
    mwIndex *R = mxGetIr(prhs[1]);
    mwIndex *starts = mxGetJc(prhs[1]);
    double* vals = mxGetPr(prhs[1]);
     
    /*Output pointers and memory*/
    double* out_color_elems = 0;
    double* out_color_starts = 0;
    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
    out_color_elems = mxGetPr(plhs[0]);

    /* The Program */
    aux_color_elems = (mwIndex*)malloc(n*sizeof(mwIndex));
    
    for ( i = 0 ; i < n ; i++){
        if (starts_t[i+1]-starts_t[i] > maximalDegree){
            maximalDegree = starts_t[i+1]-starts_t[i];
        }
    }
    if (vals!=0){
         for ( i = 0 ; i < n ; i++){
            if (starts[i+1]-starts[i] > maximalDegree){
                maximalDegree = starts[i+1]-starts[i];
            }
        }
    }
    
    aux_colors = (mwIndex*)malloc(maximalDegree*sizeof(mwIndex));
    if (aux_color_elems==0 || aux_colors == 0){
        printf("mex: multiColorSetup: Memory allocation falied!!");
        return;
    }
    for ( i = 0 ; i < maximalDegree ; i++){
        aux_colors[i] = 0;
    }
    for ( i = 0 ; i < n ; i++){
        // here aux_colors is a temporary array that marks the colors of i's neighbors
        j=starts_t[i];
        while (C_t[j]<i){
            aux_colors[aux_color_elems[C_t[j]]] = 1;
            ++j;
        }
        if (vals!=0){
            j=starts[i];
            while (R[j]<i){
                aux_colors[aux_color_elems[R[j]]] = 1;
                ++j;
            }
        }
        // At this point we know that aux_colors marks all the colors found in i's neighborhood.
        j = 0;
        while (aux_colors[j]!=0){
            aux_colors[j] = 0; 
            ++j;
        }
        aux_color_elems[i] = j;
        n_colors += (j >= n_colors);
        while (j < n_colors){
            aux_colors[j] = 0;
            ++j;
        }
	} 
    for (i = 0 ; i < n ; i++){ // count colors:
        aux_colors[aux_color_elems[i]]++;
    }
    if (n_colors > 2){
        // find two maximal && replace largest and second largest to be first (0) and last (n_colors-1)
        makeLargestTwoColorsFirstAndLast(aux_colors, aux_color_elems, n_colors, n );
    }
    plhs[1] = mxCreateDoubleMatrix(1, n_colors+1, mxREAL);
    out_color_starts = mxGetPr(plhs[1]);
    out_color_starts[0] = 0;
    // Preparing output 
    // Essentially BucketSort without the initialization:
    for (j=0; j < n_colors; ++j){
        out_color_starts[j+1] = out_color_starts[j] + aux_colors[j];
    }
    for (i=0; i < n; ++i){
        k = aux_color_elems[i];
        out_color_elems[(mwIndex)out_color_starts[k+1] - aux_colors[k]] = i;
        aux_colors[k]--;
    }
//     for (i=0, j=0; j < n_colors; ++j){
//         out_color_starts[j+1] = out_color_starts[j] + aux_colors[j];
//         for (k = aux_colors[j]; k > 0; --k){
//             out_color_elems[i++] = j;
//         }
//     }
    free(aux_colors);
    free(aux_color_elems);
}
void makeLargestTwoColorsFirstAndLast(mwIndex* aux_colors, mwIndex* aux_color_elems, mwIndex n_colors, mwIndex n ){
    mwIndex firstLargest, secondLargest, i;
    if (aux_colors[0] > aux_colors[1]){
        firstLargest = 0;
        secondLargest = 1;
    }else{
        firstLargest = 1;
        secondLargest = 0;
    }
    for (i = 2 ; i < n_colors ; i++){
        if (aux_colors[i] > aux_colors[secondLargest]){
            if (aux_colors[i] > aux_colors[firstLargest]){
                secondLargest = firstLargest ; firstLargest = i; 
            }else{
                secondLargest = i;
            }
        }
    }
    if (firstLargest!=0){
        for (i=0 ; i < n ; i++){
            if (aux_color_elems[i]==firstLargest){
                aux_color_elems[i] = 0; 
            }else if (aux_color_elems[i]==0){
                aux_color_elems[i] = firstLargest;
            }
        }
        swap(&aux_colors[0], &aux_colors[firstLargest]);
        if (secondLargest==0){
            secondLargest = firstLargest;
        }
    }
    if (secondLargest!=n_colors-1){
        for (i=0 ; i < n ; i++){
            if (aux_color_elems[i]==secondLargest){
                aux_color_elems[i] = n_colors-1; 
            }else if (aux_color_elems[i]==n_colors-1){
                aux_color_elems[i] = secondLargest;
            }
        }
        swap(&aux_colors[n_colors-1], &aux_colors[secondLargest]);
    }
//     printf("After swapping first and last:\n");
//     for (i = 0 ; i < n_colors ; i++){
//         printf("aux_colors[i] is %d\n",aux_colors[i] );
//     }
}