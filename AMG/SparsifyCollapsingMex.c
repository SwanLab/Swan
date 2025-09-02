#include "mex.h"
#include <omp.h>
#include <math.h>

void splitMatrices(mwIndex startRow, mwIndex endRow, 
        mwIndex* R_A0 , mwIndex *starts_A0, double  *V_A0 ,  
        mwIndex *R_Aomega , mwIndex *starts_Aomega, double  *V_Aomega);

void Sparsify(mwIndex startRow, mwIndex endRow, mwIndex naux,
        mwIndex* R_A0 , mwIndex *starts_A0, double  *V_A0 ,  
        mwIndex *R_Aomega , mwIndex *starts_Aomega, double *V_Aomega,
        mwIndex *R_RP0    , mwIndex *starts_RP0   , double *V_RP0,
        mwIndex *R_R0P    , mwIndex *starts_R0P   , double *V_R0P,
        mwIndex *C_A0_cols, mwIndex *starts_A0_cols,double *V_A0_cols, mwIndex* globalDiagonal);

void SparsifyAki(mwIndex k, mwIndex i, double Aki, 
        mwIndex* R_A0     , mwIndex *starts_A0     , double  *V_A0,
        mwIndex *R_RP0    , mwIndex *starts_RP0    , double *V_RP0,
        mwIndex *C_R0P    , mwIndex *starts_R0P    , double *V_R0P,
        mwIndex *C_A0_cols, mwIndex *starts_A0_cols, double *V_A0_cols,
        mwIndex *intI1I2  , mwIndex *auxI1         , mwIndex *auxI2,
        mwIndex *auxm1    , mwIndex *auxm2         , double* thetta, mwIndex* globalDiagonal ,mwIndex naux_list);


int intersect(mwIndex* A, mwIndex* B, mwIndex startA, mwIndex endA, 
        mwIndex startB, mwIndex endB, mwIndex* ans, mwIndex* iA, mwIndex* iB);
// mwIndex setDiff(mwIndex* A, mwIndex* B, mwIndex startA, mwIndex endA, 
//         mwIndex startB, mwIndex endB, mwIndex* ans, mwIndex* i);
mwIndex findBinarySearch(mwIndex* A, mwIndex startA, mwIndex endA, mwIndex toFind);
void GenerateDiagonalGlobalIndices(mwIndex startRow, mwIndex endRow, 
        mwIndex* R_A0, mwIndex* starts_A0, mwIndex* globalDiagonal);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //mex  -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims SparsifyCollapsingMex.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"
    mwIndex id, Nthrds, istart, iend;
    mwIndex naux = 0;
// // 
    mwIndex *R_A0 = mxGetIr(prhs[0]);
    mwIndex *starts_A0 = mxGetJc(prhs[0]);
    double  *V_A0 = mxGetPr(prhs[0]);
    //mwIndex nzmax_A0 = mxGetNzmax(prhs[0]);
    mwIndex n = mxGetN(prhs[0]);
    
    mwIndex *R_Aomega = mxGetIr(prhs[1]);
    mwIndex *starts_Aomega = mxGetJc(prhs[1]);
    double  *V_Aomega = mxGetPr(prhs[1]);
    mwIndex nzmax_Aomega = mxGetNzmax(prhs[1]);

    mwIndex *R_RP0 = mxGetIr(prhs[2]);
    mwIndex *starts_RP0 = mxGetJc(prhs[2]);
    double  *V_RP0 = mxGetPr(prhs[2]);

    mwIndex *C_R0P = mxGetIr(prhs[3]);
    mwIndex *starts_R0P = mxGetJc(prhs[3]);
    double  *V_R0P = mxGetPr(prhs[3]);
    
    
    mwIndex *C_A0_cols = mxGetIr(prhs[4]);
    mwIndex *starts_A0_cols = mxGetJc(prhs[4]);
    double  *V_A0_cols = mxGetPr(prhs[4]);
   
    mwIndex* globalDiagonal = (mwIndex*)malloc(n*sizeof(mwIndex));
    
    
    
    #pragma omp parallel private(id, Nthrds, istart, iend) num_threads(omp_get_num_procs())
    {
        //printf("%d",);
        id = omp_get_thread_num();
        Nthrds = omp_get_num_threads();
        istart = id * n / Nthrds;
        iend = (id+1) * n / Nthrds;
        if (id == Nthrds-1) iend = n;
        splitMatrices(istart, iend, R_A0, starts_A0, V_A0, R_Aomega, starts_Aomega, V_Aomega);
        #pragma omp barrier
                GenerateDiagonalGlobalIndices(0,n,R_A0,starts_A0,globalDiagonal);
        #pragma omp barrier
                Sparsify(istart,iend,naux, R_A0 , starts_A0, V_A0, R_Aomega, starts_Aomega, V_Aomega,
                R_RP0, starts_RP0, V_RP0, C_R0P, starts_R0P, V_R0P,
                C_A0_cols, starts_A0_cols, V_A0_cols,globalDiagonal);
		#pragma omp barrier
    }
    free(globalDiagonal);
}

void Sparsify(mwIndex startRow, mwIndex endRow, mwIndex naux,
        mwIndex* R_A0 , mwIndex *starts_A0, double  *V_A0 ,  
        mwIndex *R_Aomega , mwIndex *starts_Aomega, double *V_Aomega,
        mwIndex *R_RP0    , mwIndex *starts_RP0   , double *V_RP0,
        mwIndex *C_R0P    , mwIndex *starts_R0P   , double *V_R0P,
        mwIndex *C_A0_cols, mwIndex *starts_A0_cols,double *V_A0_cols,  mwIndex* globalDiagonal){

    mwIndex k,i, global_idx, naux_list;
    mwIndex *intI1I2, *auxI1, *auxI2, *auxm1, *auxm2;
    double *thetta;
    double Aki;
    
    for (k = startRow ; k < endRow ;++k){
        if ((starts_A0[k+1] - starts_A0[k])>naux){
            naux = (starts_A0[k+1] - starts_A0[k]);
        }
    }
	for (k = startRow ; k < endRow ;++k){
        if ((starts_A0_cols[k+1] - starts_A0_cols[k])>naux){
            naux = (starts_A0_cols[k+1] - starts_A0_cols[k]);
        }
    }
    naux = naux + 5;
	naux_list = naux*naux; 
    intI1I2 = (mwIndex*)malloc(naux*sizeof(mwIndex));
    auxI1 = (mwIndex*)malloc(naux*sizeof(mwIndex));
    auxI2 = (mwIndex*)malloc(naux*sizeof(mwIndex));
	
    auxm1 = (mwIndex*)malloc(naux_list*sizeof(mwIndex));
    auxm2 = (mwIndex*)malloc(naux_list*sizeof(mwIndex));
    thetta = (double*)malloc(naux_list*sizeof(double));
    
    for (k = startRow ; k < endRow ;++k){
        for (global_idx = starts_Aomega[k] ; global_idx < starts_Aomega[k+1]; ++global_idx){
            Aki = V_Aomega[global_idx];
            i = R_Aomega[global_idx];
            if (Aki!=0){
                //printf("Sparsifing (%d,%d,%lf)\n",k,i,Aki);
                SparsifyAki(k,i,Aki,R_A0 , starts_A0, V_A0,
                        R_RP0, starts_RP0, V_RP0, C_R0P, starts_R0P, V_R0P,
                        C_A0_cols, starts_A0_cols, V_A0_cols,
                        intI1I2, auxI1, auxI2, auxm1, auxm2, thetta, globalDiagonal,naux_list);
            }
        }
    }
    free(intI1I2); 
    free(auxI1);   
    free(auxI2);
    free(auxm1);
    free(auxm2);
    free(thetta);  
}


void SparsifyAki(mwIndex k, mwIndex i, double Aki, 
        mwIndex* R_A0     , mwIndex *starts_A0     , double  *V_A0,
        mwIndex *R_RP0    , mwIndex *starts_RP0    , double *V_RP0,
        mwIndex *C_R0P    , mwIndex *starts_R0P    , double *V_R0P,
        mwIndex *C_A0_cols, mwIndex *starts_A0_cols, double *V_A0_cols,
        mwIndex *intI1I2  , mwIndex *auxI1         , mwIndex *auxI2,
        mwIndex *auxm1    , mwIndex *auxm2         , double* thetta,  mwIndex* globalDiagonal ,mwIndex naux_list){
     
    mwIndex m1,m2, it, n_intersect, n_intersect2, m1it, list_it,
            global_km2,global_m1i,global_m1m1,global_m2m2,global_m2m1,kk;
	mwIndex *tmp_ind=0;
	double *tmp_double=0;
    double delta,R0Pm1i;
    // printf("k:%ld, i:%ld, Aki: %lf\n",k,i,Aki);
    
    // Now we check distance 2 paths between i and k
    n_intersect = intersect(C_R0P , R_RP0, starts_R0P[i] , starts_R0P[i+1],
            starts_RP0[k], starts_RP0[k+1],intI1I2, auxI1, auxI2);
    //printf("n_intersect: %d\n",n_intersect);
    if (n_intersect > 0){
        delta = 0;
        for (it = 0 ; it < n_intersect ; it++){
            thetta[it] = fabs(V_R0P[auxI1[it]]*V_RP0[auxI2[it]]);
            //printf("mex: V_R0P:%lf; V_RP0: %lf\n",V_R0P[auxI1[it]],V_RP0[auxI2[it]]);
            delta += thetta[it];
        }
        delta = Aki/delta;
        for (it = 0 ; it < n_intersect ; it++){
            thetta[it] *= delta;
        }
        for (it = 0 ; it < n_intersect ; it++){
            m1 = intI1I2[it];
            delta = thetta[it];
//             printf("mex: %ld->%ld->%ld:%lf\n",i,m1,k,delta);
            // Here: m1==m2
            // here we assume that A0 contains R0P and P0R
            global_km2 = findBinarySearch(R_A0, starts_A0[k], starts_A0[k+1], m1);
            global_m1i = findBinarySearch(R_A0, starts_A0[m1], starts_A0[m1+1], i);
            global_m1m1 = globalDiagonal[m1];
//             global_m1m1 = findBinarySearch(R_A0, starts_A0[m1], starts_A0[m1+1], m1);
            #pragma omp atomic
            V_A0[global_km2] += delta;
            #pragma omp atomic
            V_A0[global_m1i] += delta;
            #pragma omp atomic
            V_A0[global_m1m1] -= delta;
        }
    }else{
        // Now there is no distance 2 path. We seek distance 3 paths.
        list_it = 0;
        delta = 0;
        for (m1it = starts_R0P[i] ; m1it < starts_R0P[i+1] ; m1it++){
            m1 = C_R0P[m1it];
            R0Pm1i = V_R0P[m1it];
            n_intersect2 = intersect(C_A0_cols , R_RP0, starts_A0_cols[m1] , starts_A0_cols[m1+1],
                    starts_RP0[k], starts_RP0[k+1],intI1I2, auxI1, auxI2);
            for (it = 0 ; it < n_intersect2 ; it++){
                auxm1[list_it] = m1;
                auxm2[list_it] = intI1I2[it];
                thetta[list_it] = fabs(R0Pm1i*V_A0_cols[auxI1[it]]*V_RP0[auxI2[it]]);
                delta+=thetta[list_it];
                ++list_it;
				if (list_it == naux_list){
					printf("We're out of place in auxiliary arrays, but don't worry - we double naux_list\n");
					tmp_ind = (mwIndex*)malloc(2*naux_list*sizeof(mwIndex));
					for (kk=0 ; kk<naux_list ; ++kk){
						tmp_ind[kk] = auxm1[kk];
					}
					free(auxm1);
					auxm1 = tmp_ind;
					tmp_ind = (mwIndex*)malloc(2*naux_list*sizeof(mwIndex));
					for (kk=0 ; kk<naux_list ; ++kk){
						tmp_ind[kk] = auxm2[kk];
					}
					free(auxm2);
					auxm2 = tmp_ind;
					tmp_double = (double*)malloc(2*naux_list*sizeof(double));
					for (kk=0 ; kk<naux_list ; ++kk){
						tmp_double[kk] = thetta[kk];
					}
					free(thetta);
					thetta = tmp_double;
					return;
				}
            }
        }
        delta = Aki/delta;
        for (it = 0 ; it < list_it ; it++){
            thetta[it] *= delta;
        }
        for (it = 0 ; it < list_it ; it++){
            m1 = auxm1[it];
            m2 = auxm2[it];
            delta = thetta[it];
//             printf("mex: %ld->%ld->%ld->%ld:%lf\n",i,m1,m2,k,delta);
            global_km2 = findBinarySearch(R_A0, starts_A0[k], starts_A0[k+1], m2);
            global_m1i = findBinarySearch(R_A0, starts_A0[m1], starts_A0[m1+1], i);
//             global_m1m1 = findBinarySearch(R_A0, starts_A0[m1], starts_A0[m1+1], m1);
//             global_m2m2 = findBinarySearch(R_A0, starts_A0[m2], starts_A0[m2+1], m2);
            global_m2m2 = globalDiagonal[m2];
            global_m1m1 = globalDiagonal[m1];
            global_m2m1 = findBinarySearch(R_A0, starts_A0[m2], starts_A0[m2+1], m1);
            #pragma omp atomic
            V_A0[global_m2m2] -= delta;
            #pragma omp atomic
            V_A0[global_m2m1] += delta;
            #pragma omp atomic
            V_A0[global_km2]  += delta;
            #pragma omp atomic
            V_A0[global_m1i]  += delta;
            #pragma omp atomic
            V_A0[global_m1m1] -= delta;
        }
    }
}

void GenerateDiagonalGlobalIndices(mwIndex startRow, mwIndex endRow, 
        mwIndex* R_A0, mwIndex* starts_A0, mwIndex* globalDiagonal){
    mwIndex k;
    for (k = startRow ; k < endRow ;++k){
        globalDiagonal[k] = findBinarySearch(R_A0, starts_A0[k], starts_A0[k+1], k);
    }
}

void splitMatrices(mwIndex startRow, mwIndex endRow, 
        mwIndex* R_A0 , mwIndex *starts_A0, double  *V_A0 ,  
        mwIndex *R_Aomega , mwIndex *starts_Aomega, double  *V_Aomega){
    mwIndex k,global_idx0, global_idx,tmp0, tmpOmega;
    // SPLITTING
    for (k = startRow ; k < endRow ;++k){
        for (global_idx0 = starts_A0[k] ; global_idx0 < starts_A0[k+1] ; ++global_idx0){
           V_A0[global_idx0] = 0; 
        }
        global_idx0 = starts_A0[k];
        global_idx = starts_Aomega[k];
        while (global_idx < starts_Aomega[k+1] && global_idx0 < starts_A0[k+1]){
            tmp0 = R_A0[global_idx0];
            tmpOmega = R_Aomega[global_idx];
            if (tmp0 == tmpOmega){
                V_A0[global_idx0] = V_Aomega[global_idx];
                V_Aomega[global_idx] = 0;
                ++global_idx0;
                ++global_idx;
            }else{
                global_idx0 += (tmp0 < tmpOmega);
                global_idx += (tmp0 > tmpOmega);
            }
        }
    }
}

int intersect(mwIndex* A, mwIndex* B, mwIndex startA, mwIndex endA,
        mwIndex startB, mwIndex endB, mwIndex* AB, mwIndex* iA, mwIndex* iB)
{
    mwIndex ia = startA, ib = startB;
    mwIndex k = 0;
    mwIndex tmpA,tmpB;
    if (A[startA] > B[endB-1] || A[endA-1] < B[startB]){
        return k;
    }
    while ((ia < endA) && (ib < endB)){
        tmpA = A[ia];
        tmpB = B[ib];
        if (tmpA == tmpB){
            AB[k] = tmpA;
            iA[k] = ia++;
            iB[k] = ib++;
            ++k;
        }else{
            ia += tmpA < tmpB;
            ib += tmpA > tmpB;
        }
    }
    return k;
}
mwIndex findBinarySearch(mwIndex* A, mwIndex startA, mwIndex endA, mwIndex toFind){
    mwIndex i;
    for (i = startA ; i < endA ; ++i){
        if (A[i]==toFind){
           return i;
        }
    }
    printf("INDEX NOT FOUND!!! IMPOSSIBLE - THERE'S A BUG");
    return 0;
//     mwIndex mid;
//     endA = endA-1;
//     // continually narrow search until just one element remains
//     while (startA < endA)
//     {
//         mid = (startA+endA)>>1;
// //         printf("mid = %ld\n",mid);
//         // reduce the search
//         if (A[mid] < toFind){
//             startA = mid+1;
// //             printf("start = %ld\n",startA);
//         }else{
//             endA = mid;
// //             printf("end = %ld\n",endA);
//         }
//     }
//     if ((endA == startA) && (A[startA] == toFind))
//         return startA;
//     else{
//         printf("INDEX NOT FOUND!!! IMPOSSIBLE - THERE'S A BUG");
//         return 0;
//     }
}
  

// mwIndex setDiff(mwIndex* A, mwIndex* B, mwIndex startA, mwIndex endA, 
//         mwIndex startB, mwIndex endB, mwIndex* AB, mwIndex* IA){
//     mwIndex ia = startA, ib = startB;
//     mwIndex k = 0;
//     mwIndex tmpA,tmpB;
//     while ((ia < endA) && (ib < endB)){
//         tmpA = A[ia];
//         tmpB = B[ib];
//         if (tmpA < tmpB){
//             IA[k] = ++ia;
//             AB[k] = tmpA;
//             ++k;
//         }else{
//             ia += tmpA == tmpB;
//             ib += tmpA >= tmpB;
//         }
//     }
//     while (ia < endA){
//         AB[k] = A[ia];
//         IA[k] = ++ia;
//         ++k;
//     }
// }
// 