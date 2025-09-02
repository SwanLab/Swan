/*---------------------------------------------------------------------------
Copyright (2010): Eran Treister and Irad Yavneh. 
This code is distributed under the terms of the GNU General Public
License 2.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose."
---------------------------------------------------------------------------*/
#include "mex.h"
//#include<stdio.h>
#define epsilon 1e-16
class SquareMatrix{
    int N;
	int NNZ;
	int *colIdx; // usually columns
	int *rowStarts;
	double *values;
	public: 
	      int debug;
          SquareMatrix(int n, int nnz, int* starts, int* otherIdx, double* vals);
          int nnz();
          int numOfRows();
	      int getColIdx(int i);
	      double getValues(int i);
	      void setValue(int i, double v);
          int getRowStarts(int row);
	      int getRowEnds(int row);
	      int numOfNzInRow(int row);
	      void freeArrays();         
};
void multiplyDiagonalFromLeft(SquareMatrix* S,double* x);
void getMinimumsByRow(SquareMatrix* S, double* m);
void zeroLowerThanThetaMaxInRow(SquareMatrix* S, double* mins, double theta);
SquareMatrix* trimAndTranspose(SquareMatrix* S);
int nnzInSymmetrization(SquareMatrix* A, SquareMatrix* At);
void symmetrisizeAndTrim(SquareMatrix* A, SquareMatrix* At, 
        int* answerCols, double* answerValues, int* answerStarts, int N, int symNNZ);
//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int* rowStarts = (int*)mxGetData(prhs[0]);
    int* colIdx = (int*)mxGetData(prhs[1]);
    double* values = (double*)mxGetData(prhs[2]);
    int n = (int)*mxGetPr(prhs[3]);
    int nnz = (int)*mxGetPr(prhs[4]);
    double *x = (double*)mxGetData(prhs[5]);
    int *startsOut, *colsOut;
    double *valOut;
    double thetta  = (double)*mxGetPr(prhs[6]);
    SquareMatrix* S = new SquareMatrix(n,nnz, rowStarts, colIdx, values);
    multiplyDiagonalFromLeft(S,x);
    double* mins = new double[n];
    getMinimumsByRow(S, mins);
    zeroLowerThanThetaMaxInRow(S, mins, thetta);
    SquareMatrix* St = trimAndTranspose(S);
    int symNNZ = nnzInSymmetrization(S, St);
    delete mins;
    const int dims[] = {1,n};
    const int dimsNNZ[] = {1,symNNZ};
    plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    plhs[1] = mxCreateNumericArray(2,dimsNNZ,mxINT32_CLASS,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, symNNZ, mxREAL);
    startsOut = (int*)mxGetData(plhs[0]);
    colsOut = (int*)mxGetData(plhs[1]);
    valOut = mxGetPr(plhs[2]);
    symmetrisizeAndTrim(S, St, colsOut, valOut, startsOut, n, symNNZ);
    delete S;
    St->freeArrays();
    delete St;
}

//******************************************************************************
void multiplyDiagonalFromLeft(SquareMatrix* S,double* x){
    int k;
    for (k=0 ; k < S->nnz() ; k++){
        S->setValue(k,S->getValues(k)*x[S->getColIdx(k)-1]);
    }
}
void getMinimumsByRow(SquareMatrix* S, double* m){
    int k;
    double tmp;
    for (k=0 ; k < S->numOfRows() ; ++k){
        m[k] = 0.0; // here we assume that the mimimum is negative
        int rowEnds = S->getRowEnds(k+1);
        int idx = S->getRowStarts(k+1);
        while (idx <= rowEnds){
              tmp = S->getValues(idx);
              if (tmp < m[k]){
                 m[k] = tmp;
              }
              ++idx;
        }
    }
}
void zeroLowerThanThetaMaxInRow(SquareMatrix* S, double* mins, double theta){
    // this func also negates the matrix!
    int k,count = 0;
    double tmp;
    for (k=0 ; k < S->numOfRows() ; ++k){
        int rowEnds = S->getRowEnds(k+1);
        int idx = S->getRowStarts(k+1);
        while (idx <= rowEnds){
              tmp = S->getValues(idx);
              // here we assume m[k]<0
              if ( (tmp / mins[k]) < theta){
                 S->setValue(idx,0);
              }else{
                 //S->setValue(idx,-tmp);
                 // WORK WITH UNIT TEST:
                 S->setValue(idx,(tmp/mins[k]));
              }
              ++idx;
        }
    }
}

SquareMatrix* trimAndTranspose(SquareMatrix* S){
     int N = S->numOfRows();
     int* aux = new int[N];
     int idx,k, newNNZ = 0;
     int* targetStarts = new int[N];
     for (int i = 0 ; i < N ; i++){
	     aux[i] = 0; 
     }
     for (int i = 0 ; i < S->nnz() ; i++){
         if (S->getValues(i)>epsilon){
	        aux[S->getColIdx(i)-1]++;
	        newNNZ++;
         }
     }
     int* targetCols = new int[newNNZ];
     double* targetValues = new double[newNNZ];
     for (int i = 0 ; i < newNNZ ; i++){
         targetCols[i] = 0;
         targetValues[i] = 0;
     }
     idx = 0;
     for (int i = 0 ; i < N ; i++){
	     if (aux[i] == 0){
            targetStarts[i] = -1;
         }else{
            targetStarts[i] = idx;
            idx = idx + aux[i];
         }
     }
     for (int i = 1 ; i < N ; i++){
	     aux[i] += aux[i-1]; 
     }
     //aux now holds where each column ends (one after last)
     for (k=N ; k > 0 ; k--){
         int idx = S->getRowEnds(k);
         int row_k_start = S->getRowStarts(k);
         while (idx >= row_k_start){
               if (S->getValues(idx)>epsilon){  
                  int index = --aux[S->getColIdx(idx)-1];
                  targetCols[index] = k;
                  targetValues[index] = S->getValues(idx);
               }
               idx--;
         }     
     }
     SquareMatrix* transposed = new SquareMatrix(N, newNNZ, targetStarts, targetCols, targetValues);
     delete aux;
     return transposed;
     //print(transposed->values, transposed->nnz());
     //cout<<"\n";
     //printInts(transposed->colIdx, transposed->nnz());        
}

int nnzInSymmetrization(SquareMatrix* A,SquareMatrix* At){
// Count members in target symmetrized matrix 0.5*A+0.5*A':
// Assumption  - trasposed is trimmed already.
   int symNNZ = 0,k;
   int N = A->numOfRows();       
   for (k=0 ; k < N ; k++){
       int idx1 = A->getRowStarts(k+1);
       int idx2 = At->getRowStarts(k+1);
       int end1 = A->getRowEnds(k+1);
       int end2 = At->getRowEnds(k+1);
       while (idx1<=end1 && idx2<=end2){
          if (At->getColIdx(idx2) == A->getColIdx(idx1)){
              idx2++; idx1++; symNNZ++;
          }else{ 
              if (At->getColIdx(idx2) < A->getColIdx(idx1)){
                 idx2++; symNNZ++;
              }else{
                    if (A->getValues(idx1)> epsilon){
                       symNNZ++;
                    }
                    idx1++;
              }
          }        
       }
       symNNZ += (end2-idx2)+1;
       while (idx1<=end1){
           if (A->getValues(idx1)> epsilon){
                symNNZ++;
           }
           idx1++;   
       }
   }
   return symNNZ;
}

void symmetrisizeAndTrim(SquareMatrix* A, SquareMatrix* At,
                         int* answerCols, double* answerValues,
                         int* answerStarts, int N, int symNNZ){
    int idx = 0,k;
    for (k=0 ; k < N ; k++){
        int idx1 = A->getRowStarts(k+1);
        int idx2 = At->getRowStarts(k+1);
        int end1 = A->getRowEnds(k+1);
        int end2 = At->getRowEnds(k+1);
        if ((idx1!=-1) || (idx2!=-1)){
            answerStarts[k] = idx;
        }else{
            answerStarts[k] = -1;
        }
        while (idx1<=end1 && idx2<=end2){
            if (At->getColIdx(idx2) == A->getColIdx(idx1)){
                answerValues[idx] = 0.5*At->getValues(idx2) + 0.5*A->getValues(idx1);
                answerCols[idx] = A->getColIdx(idx1);
                idx2++; idx1++; idx++;
            }else{
                if (At->getColIdx(idx2) < A->getColIdx(idx1)){
                    answerValues[idx] = 0.5*At->getValues(idx2);
                    answerCols[idx] = At->getColIdx(idx2);
                    idx2++; idx++;
                }else{
                    if (A->getValues(idx1)> epsilon){
                        answerValues[idx] = 0.5*A->getValues(idx1);
                        answerCols[idx] = A->getColIdx(idx1);
                        idx++;
                    }
                    idx1++;
                }
            }
        }
        while (idx2<=end2){
            answerValues[idx] = 0.5*At->getValues(idx2);
            answerCols[idx] = At->getColIdx(idx2);
            idx2++; idx++;
        }
        while (idx1<=end1){
            if (A->getValues(idx1)> epsilon){
                answerCols[idx] = A->getColIdx(idx1);
                answerValues[idx] = 0.5*A->getValues(idx1);
                idx++;
            }
            idx1++;
        }
    }
}

//******************************************************************************
//****************** SparseMatrix **********************************************
//******************************************************************************

SquareMatrix::SquareMatrix(int n, int nnz, int* starts, int* otherIdx, double* vals) {
	N = n;
	NNZ = nnz;
	colIdx = otherIdx;
	values = vals;
	rowStarts = starts;
	debug = 1;
}
int SquareMatrix::nnz(){
	return NNZ;
}
int SquareMatrix::getColIdx(int i) {
	return colIdx[i];
}
double SquareMatrix::getValues(int i) {
	return values[i];
}
int SquareMatrix::getRowStarts(int row){
	return rowStarts[row-1];
}
int SquareMatrix::getRowEnds(int row){
	if (rowStarts[row-1]==-1)
		return -2;
	// starts[row-1+1] = starts[row]!
	while (row < N && rowStarts[row]==-1)
	{
		++row;
	}
	return (row==N) ? nnz()-1 : rowStarts[row]-1;
}
int SquareMatrix::numOfNzInRow(int row){
	return (getRowStarts(row)!=-1) ? (getRowEnds(row) - getRowStarts(row) + 1) : -1;
}
void SquareMatrix::setValue(int i, double val){
    values[i] = val;   
}
void SquareMatrix::freeArrays(){
    delete values;
    delete colIdx;
    delete rowStarts;   
}
int SquareMatrix::numOfRows(){
    return N;   
}
//******************************************************************************
//*************************     UNIT TEST     **********************************
//******************************************************************************
/*
void print(double* vals, int n){
   for (int k = 0 ; k < n ; k++){
       cout<<vals[k]<<",";
   }
   cout<<"\n";
}
void printInts(int* vals, int n){
   for (int k = 0 ; k < n ; k++){
       cout<<vals[k]<<",";
   }
   cout<<"\n";
}
int main(int argc, char *argv[])
{
       1     2     3     4     5     6     7     8     9        

1:     1   +0.33   0     0     0   -0.33 -0.25 -0.5    0
2:    -0.25  1     0   -0.5  -0.33   0     0     0     0
3:     0     0     1     0     0   +0.33 -0.25   0   -0.5
4:     0   -0.33   0     1     0     0     0   -0.5    0
5:     0   -0.33   0     0     1     0   -0.25   0   -0.5
6:    -0.25  0   -0.33   0     0     1   -0.25   0     0
7:    -0.25  0   -0.33   0   -0.33 -0.33   1     0     0
8:    +0.25  0     0   -0.5    0     0     0     1     0
9:     0     0   -0.33   0   -0.33   0     0     0     1

     int nnz = 34;
     int n = 9;
     int colIdx[] = {1,2,6,7,8,
                     1,2,4,5,
                     3,6,7,9,
                     2,4,8,
                     2,5,7,9,
                     1,3,6,7,
                     1,3,5,6,7,
                     1,4,8,
                     3,9};
     int rowStarts[] = {0,5,9,13,16,20,24,29,32};
     double vals[] = {1,+0.33,-0.33,-0.25,-0.5,
                      -0.25,1,-0.5,-0.33,
                      1,+0.33, -0.25, -0.5,
                      -0.33,1,-0.5,
                      -0.33 ,1,-0.25,-0.5,
                      -0.25,-0.33,1,-0.25,
                      -0.25 ,-0.33,-0.33, -0.33 ,  1,
                      +0.25 ,-0.5 ,1 ,
                      -0.33 ,1};
     double x[] = {1,1,1,1,1,1,1,1,1};
     int n = 5;
     int nnz = 18;
     double x[] = {1,1,1,1,1};
     int colIdx[] = {1,2,4,5,1,2,5,1,2,3,4,5,2,4,1,3,4,5};
     double vals[] = {0.65,0.6,-0.59,-0.59,-0.67,1,-0.65,-0.43,-0.14,0.25,-0.24,0.79,-0.65,0.14,-0.08,-0.97,-0.03,0.96};
     int rowStarts[] = {0,4,7,12,14};
     cout<<"Original A:";
     printInts(colIdx,nnz);
     cout<<"\n";
     print(vals,nnz);
     
     SquareMatrix* S = new SquareMatrix(n,nnz, rowStarts, colIdx, vals);     
     multiplyDiagonalFromLeft(S,x);
     double* mins = new double[n];
     getMinimumsByRow(S, mins);
     zeroLowerThanThetaMaxInRow(S, mins, 0.1);
     SquareMatrix* St = trimAndTranspose(S);
     int symNNZ = nnzInSymmetrization(S, St);
     int* answerCols = new int[symNNZ];
     double* answerValues = new double[symNNZ];
     int* answerStarts = new int[n];
     symmetrisizeAndTrim(S, St, answerCols, answerValues, answerStarts, n, symNNZ);
     
     cout<<"Symmetrized:\n";
     printInts(answerCols,symNNZ);
     cout<<"\n";
     print(answerValues,symNNZ);
     delete S;
     delete mins;
     St->freeArrays();
     delete St;
     
     delete answerCols;
     delete answerValues;
     delete answerStarts;

     system("PAUSE");
     return EXIT_SUCCESS;
}
*/
