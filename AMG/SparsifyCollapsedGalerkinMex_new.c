#include "mex.h"
#include <omp.h>
#include <math.h>

//cd .\MEXfunc ; mex -O 'CXXOPTIMFLAGS=-DNDEBUG -O3' -largeArrayDims SparsifyCollapsedGalerkinMex_new.c COMPFLAGS="$COMPFLAGS -openmp" LINKFALGS="$LINKFALGS -openmp"; cd .\..
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
typedef struct member_t{
    double v_vec;
    mwIndex i_vec;
    mwIndex j_mat;    
    mwIndex global_mat;
} member;

int FirstSmaller(member* A, member* B){
    return (A->j_mat < B->j_mat);
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
typedef struct Heap_t{
    int heapSize;
    int maxSize;
    member** arr;
} Heap;

void Insert(Heap* heap, member* element)
{
    int par;
    int now = heap->heapSize;
    if (now == heap->maxSize-1){
		
        printf("ERROR: Heap not big enough for size %d Increase MAX_HEAP_SIZE)!\n",heap->maxSize);
        return;
    }
    heap->heapSize++;
    while (now > 0){
        par = (now - 1)>>1;
        if (FirstSmaller(element,heap->arr[par])){
            heap->arr[now] = heap->arr[par];
            now = par;
        }else{
            break;
        }
    }
    heap->arr[now] = element;
}

member* DeleteMin(Heap* heap)
{
    int l,r,child,k=0;
    member* least = heap->arr[k];
    member* x;
    if (heap->heapSize==0){
        return NULL;
    }
    
    --heap->heapSize;
    x = heap->arr[heap->heapSize];
    while (1) {
        l = 2 * k + 1;;
        if (l >= heap->heapSize) {
            break;
        }else {
            r = 2 * (k + 1);
            child = (r >= heap->heapSize || (FirstSmaller(heap->arr[l],heap->arr[r]))) ? l : r;
            if (FirstSmaller(heap->arr[child],x)) {
                heap->arr[k] = heap->arr[child];
                k = child;
            }
            else{
                break;
            }
        }
    }
    heap->arr[k] = x;
    //  Prevent leakage...?
    return least;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

int intersect(mwIndex* A, mwIndex* B, mwIndex startA, mwIndex endA,
        mwIndex startB, mwIndex endB, mwIndex* AB, mwIndex* iA, mwIndex* iB);
mwIndex findLinearSearch(mwIndex* A, mwIndex startA, mwIndex endA, mwIndex toFind);
mwIndex findBinarySearch(mwIndex* A, mwIndex startA, mwIndex endA, mwIndex toFind);


mwIndex MultyplyRowTimesMat(mwIndex *vec_ind , double *vec_vals, mwIndex nnz_vec,
        mwIndex *C_mat  , mwIndex *starts_mat , double *V_mat,
        mwIndex *ans_ind , double *ans_vals, Heap* heap, member* members);

int SparsifyAki(mwIndex k, mwIndex i, double Aki, 
        mwIndex* R_A0     , mwIndex *starts_A0     , double  *V_A0,
        mwIndex *R_RP0    , mwIndex *starts_RP0    , double *V_RP0,
        mwIndex *C_R0P    , mwIndex *starts_R0P    , double *V_R0P,
        mwIndex *C_A0_cols, mwIndex *starts_A0_cols, double *V_A0_cols,
        mwIndex *intI1I2  , mwIndex *auxI1         , mwIndex *auxI2,
        mwIndex *auxm1    , mwIndex *auxm2         , double* thetta,  mwIndex* globalDiagonal , int MAX_HEAP_SIZE_nc);


void MultyplyAndSparsify(mwIndex startRow, mwIndex endRow, mwIndex n,
        mwIndex* C_A0 , mwIndex *starts_A0, double *V_A0,
        mwIndex *C_R  , mwIndex *starts_R , double *V_R,
        mwIndex *C_A  , mwIndex *starts_A , double *V_A,
        mwIndex *C_P  , mwIndex *starts_P , double *V_P,
        mwIndex *R_RP0    , mwIndex *starts_RP0   , double *V_RP0,
        mwIndex *C_R0P    , mwIndex *starts_R0P   , double *V_R0P,
        mwIndex *C_A0_cols, mwIndex *starts_A0_cols,double *V_A0_cols,
        mwIndex* globalDiagonal, int MAX_HEAP_SIZE_nc,int MAX_HEAP_SIZE_n){
    
    mwIndex k,global_idx,nnz_in_row,i, *aux_ind , global_idx0,tmp0, tmpOmega , sparsify_indcator,ans;
    mwIndex *intI1I2, *auxI1, *auxI2, *auxm1, *auxm2;
    double *thetta;
    double Aki;
    
    
    member* members;
    double *aux_vals;
    Heap heap;
    heap.heapSize = 0;
    heap.maxSize = MAX_HEAP_SIZE_n*MAX_HEAP_SIZE_n*MAX_HEAP_SIZE_n;
    
    
    members = (member*)malloc(heap.maxSize*sizeof(member));
    heap.arr = (member**)malloc(heap.maxSize*sizeof(member*)); 
    aux_ind = (mwIndex*)malloc(MAX_HEAP_SIZE_n*MAX_HEAP_SIZE_n*MAX_HEAP_SIZE_n*sizeof(mwIndex));
    aux_vals = (double*)malloc(MAX_HEAP_SIZE_n*MAX_HEAP_SIZE_n*MAX_HEAP_SIZE_n*sizeof(double));
   
    intI1I2 = (mwIndex*)malloc(MAX_HEAP_SIZE_nc*sizeof(mwIndex));
    auxI1 = (mwIndex*)malloc(MAX_HEAP_SIZE_nc*sizeof(mwIndex));
    auxI2 = (mwIndex*)malloc(MAX_HEAP_SIZE_nc*sizeof(mwIndex));
    auxm1 = (mwIndex*)malloc(MAX_HEAP_SIZE_nc*MAX_HEAP_SIZE_nc*MAX_HEAP_SIZE_nc*sizeof(mwIndex));
    auxm2 = (mwIndex*)malloc(MAX_HEAP_SIZE_nc*MAX_HEAP_SIZE_nc*MAX_HEAP_SIZE_nc*sizeof(mwIndex));
    thetta = (double*)malloc(MAX_HEAP_SIZE_nc*MAX_HEAP_SIZE_nc*MAX_HEAP_SIZE_nc*sizeof(double));
    
    
    
    for (k = startRow ; k < endRow ;++k){
        global_idx = starts_R[k];
        nnz_in_row = starts_R[k+1] - starts_R[k];
		if (nnz_in_row > heap.maxSize){
			printf("heap size is set to %d and we need %d",heap.maxSize,nnz_in_row);
		}
        nnz_in_row = MultyplyRowTimesMat(C_R+global_idx, V_R+global_idx, nnz_in_row,
        C_A, starts_A, V_A, aux_ind , aux_vals, &heap, members);
		if (nnz_in_row > heap.maxSize){
			printf("heap size is set to %d and we need %d",heap.maxSize,nnz_in_row);
		}
        nnz_in_row = MultyplyRowTimesMat(aux_ind, aux_vals, nnz_in_row,
        C_P, starts_P, V_P, aux_ind , aux_vals, &heap, members);
           
//         for (i = 0 ; i < nnz_in_row ;++i){
//             V_A0[starts_A0[k] + i] = aux_vals[i];
//             C_A0[starts_A0[k] + i] = aux_ind[i];
//         }
        // SPLITTING
        
        global_idx0 = starts_A0[k];
        i = 0;
        while (i<nnz_in_row && global_idx0 < starts_A0[k+1]){
            tmp0 = C_A0[global_idx0];
            tmpOmega = aux_ind[i];
            if (tmp0 == tmpOmega){
                 #pragma omp atomic
                V_A0[global_idx0] += aux_vals[i];
                global_idx0++;
                i++;
//                 printf("exit 1\n");
            }else{
                if (tmp0 < tmpOmega){
                    global_idx0++;
//                     printf("exit 2\n");
                }else{
					ans = SparsifyAki(k,aux_ind[i],aux_vals[i],C_A0 , starts_A0, V_A0,
                        R_RP0, starts_RP0, V_RP0, C_R0P, starts_R0P, V_R0P,
                        C_A0_cols, starts_A0_cols, V_A0_cols,
                        intI1I2, auxI1, auxI2, auxm1, auxm2, thetta, globalDiagonal,MAX_HEAP_SIZE_nc);
//                     printf("Sparsifing (%d,%d,%lf)\n",k,aux_ind[i],aux_vals[i]);
                    i++;
                }
            }
        }
        while (i<nnz_in_row){
			ans = SparsifyAki(k,aux_ind[i],aux_vals[i],C_A0 , starts_A0, V_A0,
                    R_RP0, starts_RP0, V_RP0, C_R0P, starts_R0P, V_R0P,
                    C_A0_cols, starts_A0_cols, V_A0_cols,
                    intI1I2, auxI1, auxI2, auxm1, auxm2, thetta, globalDiagonal,MAX_HEAP_SIZE_nc); 
//             printf("Sparsifing (%d,%d,%lf)\n",k,aux_ind[i],aux_vals[i]);
            i++;
        }
        
    }

    free(aux_ind);
    free(aux_vals);
    free(members);
    free(heap.arr);
    
    free(intI1I2); 
    free(auxI1);   
    free(auxI2);
    free(auxm1);
    free(auxm2);
    free(thetta);  
}

int SparsifyAki(mwIndex k, mwIndex i, double Aki, 
        mwIndex* R_A0     , mwIndex *starts_A0     , double  *V_A0,
        mwIndex *R_RP0    , mwIndex *starts_RP0    , double *V_RP0,
        mwIndex *C_R0P    , mwIndex *starts_R0P    , double *V_R0P,
        mwIndex *C_A0_cols, mwIndex *starts_A0_cols, double *V_A0_cols,
        mwIndex *intI1I2  , mwIndex *auxI1         , mwIndex *auxI2,
        mwIndex *auxm1    , mwIndex *auxm2         , double* thetta,  mwIndex* globalDiagonal ,mwIndex MAX_HEAP_SIZE_nc){
     
    mwIndex m1,m2, it, n_intersect, n_intersect2, m1it, list_it,
            global_km2,global_m1i,global_m1m1,global_m2m2,global_m2m1;
    double delta,R0Pm1i;
    int ans = 0;
    // printf("k:%ld, i:%ld, Aki: %lf\n",k,i,Aki);
    if (fabs(Aki)<1e-14){
		ans = 1;
		return;
	}
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
            ans = 1;
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
                delta += thetta[list_it];
                ++list_it;
				if (list_it == MAX_HEAP_SIZE_nc){
					printf("Error: we're out of place in auxiliary arrays, increase MAX_HEAP_SIZE_nc\n");
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
            //printf("mex: %ld->%ld->%ld->%ld:%lf\n",i,m1,m2,k,delta);
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
            ans = 1;
        }
    }
	if (ans==0){
        printf("Element (%ld,%ld,%.16lf) not sparsified!!!\n",k,i,Aki);
	}
}

mwIndex MultyplyRowTimesMat(mwIndex *vec_ind , double *vec_vals, mwIndex nnz_vec,
        mwIndex *C_mat  , mwIndex *starts_mat , double *V_mat,
        mwIndex *ans_ind , double *ans_vals, Heap* heap, member* members){
    mwIndex k;
    double ans_tmp_v = 0.0;
    
    mwIndex ans_tmp_global = 0;
    member *curr_min;
    heap->heapSize = 0;
    if (nnz_vec == 0){
        return ans_tmp_global;
    }
    for (k = 0 ; k < nnz_vec ; ++k){
        curr_min = &members[k];
        curr_min->v_vec = vec_vals[k];
        curr_min->i_vec = vec_ind[k];
        curr_min->global_mat = starts_mat[vec_ind[k]];
        curr_min->j_mat = C_mat[curr_min->global_mat];
        Insert(heap,curr_min);
    }
    // here we're done reading from vec_ind. It will be run over if vec_ind==ans_ind
    curr_min = DeleteMin(heap);
    ans_ind[0] = curr_min->j_mat;
    while (curr_min!=NULL){
        if (curr_min->j_mat == ans_ind[ans_tmp_global]){
//             // here we add to the answer
            ans_tmp_v += curr_min->v_vec*V_mat[curr_min->global_mat];
        }else{
            ans_vals[ans_tmp_global] = ans_tmp_v;
            ans_tmp_global++;
            ans_tmp_v = curr_min->v_vec*V_mat[curr_min->global_mat];
            ans_ind[ans_tmp_global] = curr_min->j_mat;
        }
        curr_min->global_mat++;
        if (curr_min->global_mat < starts_mat[curr_min->i_vec+1]){
            curr_min->j_mat = C_mat[curr_min->global_mat];
            Insert(heap,curr_min);
        }
        curr_min = DeleteMin(heap);
    }
    ans_vals[ans_tmp_global]=ans_tmp_v;
    ans_tmp_global++;
    return ans_tmp_global;
}

void GenerateDiagonalGlobalIndices(mwIndex startRow, mwIndex endRow, 
        mwIndex* C_A0, mwIndex* starts_A0, mwIndex* globalDiagonal){
    mwIndex k;
    for (k = startRow ; k < endRow ;++k){
        globalDiagonal[k] = findLinearSearch(C_A0, starts_A0[k], starts_A0[k+1], k);
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    mwIndex id, Nthrds, istart, iend,k;
    
    mwIndex *C_A0 = mxGetIr(prhs[0]);
    mwIndex *starts_A0 = mxGetJc(prhs[0]);
    double  *V_A0 = mxGetPr(prhs[0]);
    mwIndex nzmax_A0 = mxGetNzmax(prhs[0]);
    mwIndex nc = mxGetN(prhs[0]);
      
    mwIndex *C_R = mxGetIr(prhs[1]);
    mwIndex *starts_R = mxGetJc(prhs[1]);
    double  *V_R = mxGetPr(prhs[1]);
    
    mwIndex *C_A = mxGetIr(prhs[2]);
    mwIndex *starts_A = mxGetJc(prhs[2]);
    double  *V_A = mxGetPr(prhs[2]);
	mwIndex n = mxGetN(prhs[2]);
    
    mwIndex *C_P = mxGetIr(prhs[3]);
    mwIndex *starts_P = mxGetJc(prhs[3]);
    double  *V_P = mxGetPr(prhs[3]);
    
    mwIndex *R_RP0 = mxGetIr(prhs[4]);
    mwIndex *starts_RP0 = mxGetJc(prhs[4]);
    double  *V_RP0 = mxGetPr(prhs[4]);

    mwIndex *C_R0P = mxGetIr(prhs[5]);
    mwIndex *starts_R0P = mxGetJc(prhs[5]);
    double  *V_R0P = mxGetPr(prhs[5]);
    
    
    mwIndex *C_A0_cols = mxGetIr(prhs[6]);
    mwIndex *starts_A0_cols = mxGetJc(prhs[6]);
    double  *V_A0_cols = mxGetPr(prhs[6]);

    mwIndex* globalDiagonal = (mwIndex*)malloc(nc*sizeof(mwIndex));
	
	mwIndex MAX_HEAP_SIZE_nc = 0;
	mwIndex MAX_HEAP_SIZE_n = 0;
	
	for (k = 0 ; k < nc ;++k){
        if ((starts_A0[k+1] - starts_A0[k]) > MAX_HEAP_SIZE_nc){
            MAX_HEAP_SIZE_nc = (starts_A0[k+1] - starts_A0[k]);
        }   
    }
	for (k = 0 ; k < n ;++k){
        if ((starts_A[k+1] - starts_A[k])>MAX_HEAP_SIZE_n){
            MAX_HEAP_SIZE_n = (starts_A[k+1] - starts_A[k]);
        }   
    }
	
	
    
    #pragma omp parallel private(id, Nthrds, istart, iend, k) num_threads(omp_get_num_procs()/2)
    {
        id = omp_get_thread_num();
        //printf("num threads: %u, ",omp_get_num_threads());
        Nthrds = omp_get_num_threads();
        istart = id * nzmax_A0 / Nthrds;
        iend = (id+1) * nzmax_A0 / Nthrds;
        if (id == Nthrds-1) iend = nzmax_A0;
        for (k = istart ; k < iend ; ++k){
            V_A0[k] = 0.0;
        }
        #pragma omp barrier
        istart = id * nc / Nthrds;
        iend = (id+1) * nc / Nthrds;
        if (id == Nthrds-1) iend = nc;
        GenerateDiagonalGlobalIndices(istart,iend,C_A0,starts_A0,globalDiagonal);
        #pragma omp barrier
        MultyplyAndSparsify(istart, iend, nc, C_A0 , starts_A0, V_A0,
        C_R, starts_R, V_R, C_A  , starts_A , V_A, C_P, starts_P, V_P, 
        R_RP0, starts_RP0, V_RP0, C_R0P, starts_R0P, V_R0P, C_A0_cols, starts_A0_cols, V_A0_cols,
        globalDiagonal, MAX_HEAP_SIZE_nc,MAX_HEAP_SIZE_n);
		#pragma omp barrier		
    }   
    free(globalDiagonal);
}



mwIndex findLinearSearch(mwIndex* A, mwIndex startA, mwIndex endA, mwIndex toFind){
    mwIndex mid;
    for (mid = startA ; mid < endA ; ++mid){
        if (A[mid]==toFind){
            return mid;
        }
    }
    printf("LINEAR: INDEX NOT FOUND!!! IMPOSSIBLE - THERE'S A BUG see %d\n",toFind);
    return 0;
}


mwIndex findBinarySearch(mwIndex* A, mwIndex startA, mwIndex endA, mwIndex toFind){
    mwIndex mid;
    if (endA - startA < 50){
		mid = findLinearSearch( A, startA, endA, toFind);
		return mid;
    }else{
        endA--;
        // continually narrow search until just one element remains
        while (startA < endA)
        {
            mid = (startA+endA)>>1;
            // reduce the search
            if (A[mid] < toFind){
                startA = mid+1;
            }else{
                endA = mid;
            }
        }
        if (A[startA] == toFind)
            return startA;
        else{
            printf("BINARY: INDEX NOT FOUND!!! IMPOSSIBLE - THERE'S A BUG\n");
            return 0;
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
