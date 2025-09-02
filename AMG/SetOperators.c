#include "mex.h"
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int ia=0,ib=0,k=1;
    double *n_ans;
    int nA = max(mxGetM(prhs[0]),mxGetN(prhs[0]));
    double *A = mxGetPr(prhs[0]);
    int nB = max(mxGetM(prhs[1]),mxGetN(prhs[1]));
    double *B = mxGetPr(prhs[1]);
    double *IA,*IB,*AB;
    double tmpA,tmpB;
    char op = *((char*)mxGetData(prhs[2]));
    // we assume that arrays A and B are sorted from small to large
    if (op=='u'){ // union
        plhs[1] = mxCreateDoubleMatrix(1,nA, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(1,nB, mxREAL);
        IA = mxGetPr(plhs[1]);
        IB = mxGetPr(plhs[2]);
        /* Assign pointers to each input and output. */
        while ((ia < nA) && (ib < nB)){
            tmpA = A[ia];
            tmpB = B[ib];
            IA[ia] = k;
            IB[ib] = k;
            ia += tmpA<=tmpB;
            ib += tmpA>=tmpB;
            ++k;
        }
        while (ia<nA){
            //AB(k-1) = A(ia);
            IA[ia] = k;
            ++ia;
            ++k;
        }
        while (ib<nB){
            //AB(k-1) = B(ib);
            IB[ib] = k;
            ++ib;
            ++k;
        }
        plhs[0] = mxCreateDoubleMatrix(1,k-1, mxREAL);
        AB = mxGetPr(plhs[0]);
        for (ia = 0 ; ia < nA ; ++ia){
            AB[(int)IA[ia]-1] = A[ia];
        }
        for (ib = 0 ; ib < nB ; ++ib){
            AB[(int)IB[ib]-1] = B[ib];
        }
    }
    if (op=='i'){ // intersect
        plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
        k = 0;
        if (A[0] > B[nB-1] || A[nA-1] < B[0]){
            plhs[0] = mxCreateDoubleMatrix(1,0, mxREAL);
            plhs[1] = mxCreateDoubleMatrix(1,0, mxREAL);
            plhs[2] = mxCreateDoubleMatrix(1,0, mxREAL);
            return;
        }
        k = 0; ia = 0; ib = 0;
        while ((ia < nA) && (ib < nB)){
            tmpA = A[ia];
            tmpB = B[ib];
            k  += tmpA == tmpB;
            ia += tmpA <= tmpB;
            ib += tmpA >= tmpB;
        }
//         k = min(nB,nA);
        plhs[0] = mxCreateDoubleMatrix(1,k, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(1,k, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(1,k, mxREAL);
        AB = mxGetPr(plhs[0]);
        IA = mxGetPr(plhs[1]);
        IB = mxGetPr(plhs[2]);
        n_ans = mxGetPr(plhs[3]);
        k = 0; ia = 0; ib = 0;
        while ((ia < nA) && (ib < nB)){
            tmpA = A[ia];
            tmpB = B[ib];
            if (tmpA == tmpB){
                AB[k] = tmpA;
                IA[k] = ++ia;
                IB[k] = ++ib;
                k++;
            }else{
                ia += tmpA < tmpB;
                ib += tmpA > tmpB;
            }
            
        }
        *n_ans = k;
        return;
    }

    if (op=='d'){
        k = 0; ia = 0; ib = 0;
        while ((ia < nA) && (ib < nB)){
            tmpA = A[ia];
            tmpB = B[ib];
            if (tmpA < tmpB){
                ++k;
                ++ia;
            }else{
                ia += tmpA == tmpB;
                ib += tmpA >= tmpB;
            }
        }
        k += nA-ia;
        plhs[0] = mxCreateDoubleMatrix(1,k, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(1,k, mxREAL);
        AB = mxGetPr(plhs[0]);
        IA = mxGetPr(plhs[1]);
        k = 0; ia = 0; ib = 0;
        while ((ia < nA) && (ib < nB)){
            tmpA = A[ia];
            tmpB = B[ib];
            if (tmpA < tmpB){
                IA[k] = ++ia;
                AB[k] = tmpA;
                ++k;
            }else{
                ia += tmpA == tmpB;
                ib += tmpA >= tmpB;
            }
        }
        while (ia<nA){
            AB[k] = A[ia];
            IA[k] = ++ia;
            ++k;
        }
    }
}