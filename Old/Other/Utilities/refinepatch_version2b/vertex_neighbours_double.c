#include "mex.h"
#include "math.h"

/*
 * Ne=vertex_neighbours_double2(Fa,Fb,Fc,Vx,Vy,Vz)
 *
 */

/* Coordinates to index */
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}
int mindex2(int x, int y, int sizx) { return y*sizx+x;}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* All inputs */
    double *FacesA, *FacesB, *FacesC, *VerticesX, *VerticesY, *VerticesZ;
    /* Unsorted Neighborh list */
    int **NU, *NU_length;
    
    /* Neighbour cell array (output) */
    mxArray *Ne;
    
    
    /* Neighbour sort list of one vertex */
    int *Pneig, PneighPos, PneighStart, *Pneighf;
    double *Pneigd;
    mxArray *Pneig_matlab;
    
    /* Number of faces */
    const mwSize *FacesDims;
    int FacesN=0;
    
    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN=0;
    
    
    /* Loop variable */
    int i, j, index1, index2;
    
    /* Found */
    int found, found2;
    
    /* face vertices int */
    int vertexa, vertexb, vertexc;
    
    /* neighbour cell array length (same as vertices) */
    mwSize outputdims[1]={0};
    
    /* Check for proper number of arguments. */
    if(nrhs!=6) {
        mexErrMsgTxt("6 inputs are required.");
    } else if(nlhs!=1) {
        mexErrMsgTxt("1 output is required");
    }
    
    
    /* Read all inputs (faces and vertices) */
    FacesA=mxGetPr(prhs[0]);
    FacesB=mxGetPr(prhs[1]);
    FacesC=mxGetPr(prhs[2]);
    VerticesX=mxGetPr(prhs[3]);
    VerticesY=mxGetPr(prhs[4]);
    VerticesZ=mxGetPr(prhs[5]);
    
    /* Get number of FacesN */
    FacesDims = mxGetDimensions(prhs[0]);
    FacesN=FacesDims[0]*FacesDims[1];
    
    /* Get number of VertexN */
    VertexDims = mxGetDimensions(prhs[3]);
    VertexN=VertexDims[0]*VertexDims[1];
    
    outputdims[0]=VertexN;
    plhs[0]=mxCreateCellArray(1, outputdims);
    Ne=plhs[0];
    
    /* Neighborhs Unsorted */
    NU = (int **)malloc( VertexN* sizeof(int*) );
    NU_length = (int *)malloc( VertexN* sizeof(int) );
    for (i=0; i<VertexN; i++) {
        /* Set number of vertex neighbors to zero */
        NU_length[i]=0;
        /* Make room for neighbors */
        NU[i]=(int *)malloc( 20 * sizeof(int ) );
    }
    
    
    /* Loop throuh all faces */
    for (i=0; i<FacesN; i++) {
        /* Add the neighbors of each vertice of a face
         * to his neighbors list. */
        
        vertexa=(int)FacesA[i]-1; vertexb=(int)FacesB[i]-1; vertexc=(int)FacesC[i]-1;
        
        if(NU_length[vertexa]>18) { /* If need make extra room to store values */
            NU[vertexa]=(int *) realloc( NU[vertexa], (NU_length[vertexa]+2) * sizeof(int));
        }
        NU[vertexa][NU_length[vertexa]]  = vertexb+1; NU[vertexa][NU_length[vertexa]+1]= vertexc+1;
        NU_length[vertexa]+=2;
        
        if(NU_length[vertexb]>18) { /* If need make extra  store values */
            NU[vertexb]=(int *) realloc( NU[vertexb], (NU_length[vertexb]+2) * sizeof(int));
        }
        NU[vertexb][NU_length[vertexb]]  = vertexc+1; NU[vertexb][NU_length[vertexb]+1]=  vertexa+1;
        NU_length[vertexb]+=2;
        
        if(NU_length[vertexc]>18) { /* If need make extra store values */
            NU[vertexc]=(int *) realloc( NU[vertexc], (NU_length[vertexc]+2) * sizeof(int));
        }
        NU[vertexc][NU_length[vertexc]]  =  vertexa+1; NU[vertexc][NU_length[vertexc]+1]=  vertexb+1;
        NU_length[vertexc]+=2;
    }
    
    /*  Loop through all neighbor arrays and sort them (Rotation same as faces) */
    
    for (i=0; i<VertexN; i++) {
        if(NU_length[i]>0) {
            
            /* Create Matlab array for sorted neighbours of vertex*/
            Pneig=(int *)malloc(NU_length[i]*sizeof(int));
            PneighPos=0;
            Pneighf=NU[i];
            
            /* Start with the first vertex or if exist with a unique vertex */
            PneighStart=0;
            for(index1=0; index1<NU_length[i]; index1+=2) {
                found=0;
                for(index2=1; index2<NU_length[i]; index2+=2) {
                    if(Pneighf[index1]==Pneighf[index2]) {
                        found=1; break;
                    }
                }
                if(found==0) {
                    PneighStart=index1; break;
                }
            }
            
            Pneig[PneighPos]=Pneighf[PneighStart];    PneighPos++;
            Pneig[PneighPos]=Pneighf[PneighStart+1]; PneighPos++;
            
            
            /* Add the neighbours with respect to original rotation */
            for(j=1+found; j<(NU_length[i]/2); j++) {
                found=0;
                for(index1=0; index1<NU_length[i]; index1+=2) {
                    if(Pneighf[index1]==Pneig[PneighPos-1]) {
                        found2=0;
                        for(index2=0; index2<PneighPos; index2++) {
                            if(Pneighf[index1+1]==Pneig[index2]) { found2=1; }
                        }
                        if(found2==0) {
                            found=1;
                            Pneig[PneighPos]=Pneighf[index1+1];  PneighPos++;
                        }
                    }
                }
                if(found==0) /* This only happens with weird edge vertices */
                {
                    for(index1=0; index1<NU_length[i]; index1+=2) {
                        found2=0;
                        for(index2=0; index2<PneighPos; index2++) {
                            if(Pneighf[index1]==Pneig[index2]) { found2=1; }
                        }
                        if(found2==0) {
                            Pneig[PneighPos]=Pneighf[index1];  PneighPos++;
                            if(Pneighf[index1]==Pneig[PneighPos-1]) {
                                found2=0;
                                for(index2=0; index2<PneighPos; index2++) {
                                    if(Pneighf[index1+1]==Pneig[index2]) { found2=1; }
                                }
                                if(found2==0) {
                                    found=1;
                                    Pneig[PneighPos]=Pneighf[index1+1];  PneighPos++;
                                }
                            }
                        }
                        
                    }
                }
            }
            
            /* Add forgotten neigbours */
            if(PneighPos<NU_length[i]) {
                for(index1=0; index1<NU_length[i]; index1++) {
                    found2=0;
                    for(index2=0; index2<PneighPos; index2++) {
                        if(Pneighf[index1]==Pneig[index2]) { found2=1; break;}
                    }
                    if(found2==0) {
                        Pneig[PneighPos]=Pneighf[index1];  PneighPos++;
                    }
                }
            }
            
            outputdims[0]=PneighPos; /*(NU_length[i]/2);  */
            
            Pneig_matlab=mxCreateNumericArray(1, outputdims, mxDOUBLE_CLASS, mxREAL);
            Pneigd=(double*)mxGetPr(Pneig_matlab);

            
            /* Copy int to double array */
            for(j=0; j<outputdims[0]; j++) { Pneigd[j]=(double)Pneig[j]; }
            free(Pneig);
            mxSetCell(Ne, i, mxDuplicateArray(Pneig_matlab));
        }
    }
    
    /* Free memory */
    for (i=0; i<VertexN; i++) { free(NU[i]); }
    free(NU);
    
}

