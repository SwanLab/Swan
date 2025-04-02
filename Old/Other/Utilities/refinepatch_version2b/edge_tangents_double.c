#include "mex.h"
#include "math.h"

/*
 * function [ET_table,EV_table,ETV_index]=edge_tangents(V,Ne)
 */

/* Coordinates to index */
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}
int mindex2(int x, int y, int sizx) { return y*sizx+x;}

__inline double pow2(double val){ return val*val; }


/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* Vertex list input */
    double *V;
    /* Neighbour list input */
    const mxArray *Ne;
    mxArray *PneigMatlab;
    mwSize *PneigDims;
    int PneigLenght=0;
    
    /* Outputs */
    double *ET_table; /* Edge tangents table */
    double *EV_table; /* Edge velocity */
    double *ETV_index; /* Edge tangents/velocity index for tables */
    
    double Ea, Eb, Ec, s, h, x, Tb3D_length, Vv;
    double Np[3], Ns[3], Tb[3], Pm[3], X3[3], Y3[3], Tb3D[3];
    double *Pn, *Pnop, *Pneig;
    /* Current vertex */
    double P[3]={0,0,0};
    
    /* Table Size */
    mwSize table_Dims[2]={1,3};
    
    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN=0;
    
    /* Number of indices */
    int ETV_num=0;
    
    /* Neighbour indices */
    double neg;
    int neg1, neg2;
    
    /* Loop variables */
    int i, j, k;
    
    /* Check for proper number of arguments. */
    if(nrhs!=2) {
        mexErrMsgTxt("2 inputs are required.");
    } else if(nlhs!=3) {
        mexErrMsgTxt("3 output is required");
    }
    
    
    /* Connect Inputs */
    V=(double *)mxGetPr(prhs[0]);
    Ne=prhs[1];
    
    /* Get number of VertexN */
    VertexDims = mxGetDimensions(prhs[0]);
    VertexN=VertexDims[0];
    
    
    /* Reserve memory */
    table_Dims[0]=VertexN*6; table_Dims[1]=3;
    plhs[0]= mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    table_Dims[0]=VertexN*6; table_Dims[1]=1;
    plhs[1]= mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    table_Dims[0]=VertexN*6; table_Dims[1]=2;
    plhs[2] = mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    
    /* Connect Outputs */
    ET_table=(double *)mxGetPr(plhs[0]);
    EV_table=(double *)mxGetPr(plhs[1]);
    ETV_index=(double *)mxGetPr(plhs[2]);
    
    
    
    Pn= (double *)malloc(VertexN*3*sizeof(double));
    Pnop=(double *) malloc(VertexN*3*sizeof(double));
    

    for (i=0; i<VertexN; i++) {
        P[0]=V[i]; P[1]=V[i+VertexN]; P[2]=V[i+VertexN*2];
        PneigMatlab=mxGetCell(Ne, i);
        if( PneigMatlab == NULL)
        {
            PneigLenght=-1;
        }
        else
        {
            PneigDims=(mwSize *)mxGetDimensions(PneigMatlab);
            PneigLenght=(PneigDims[0]*PneigDims[1]);
            Pneig=(double *)mxGetPr(PneigMatlab);
        }
        /* Find the opposite vertex of each neigbourh vertex. */
        /* incase of odd number of neigbourhs interpolate the opposite neigbourh */
        if((PneigLenght%2)==0) {
            for (k=0; k<PneigLenght; k++) {
                neg =k+((double)PneigLenght)/2;
                neg1=(int)floor(neg);
                if(neg1>(PneigLenght-1)) { neg1=neg1-PneigLenght; }
                Pn[k]          =V[(int)Pneig[k]-1];
                Pn[k+VertexN]  =V[(int)Pneig[k]-1+VertexN];
                Pn[k+VertexN*2]=V[(int)Pneig[k]-1+VertexN*2];
                
                Pnop[k]          =V[(int)Pneig[neg1]-1];
                Pnop[k+VertexN]  =V[(int)Pneig[neg1]-1+VertexN];
                Pnop[k+VertexN*2]=V[(int)Pneig[neg1]-1+VertexN*2];
            }
        }
        else {
            for (k=0; k<PneigLenght; k++) {
                neg=(double)k+((double)PneigLenght)*0.5;
                neg1=(int)floor(neg); neg2=(int)ceil(neg);
                if(neg1>(PneigLenght-1)) { neg1=neg1-PneigLenght; }
                if(neg2>(PneigLenght-1)) { neg2=neg2-PneigLenght; }
                Pn[k]          =V[(int)Pneig[k]-1];
                Pn[k+VertexN]  =V[(int)Pneig[k]-1+VertexN];
                Pn[k+VertexN*2]=V[(int)Pneig[k]-1+VertexN*2];
                
                Pnop[k]          =(V[(int)Pneig[neg1]-1]           + V[(int)Pneig[neg2]-1])/2;
                Pnop[k+VertexN]  =(V[(int)Pneig[neg1]-1+VertexN]   + V[(int)Pneig[neg2]-1+VertexN])/2;
                Pnop[k+VertexN*2]=(V[(int)Pneig[neg1]-1+VertexN*2] + V[(int)Pneig[neg2]-1+VertexN*2])/2;
            }
        }
        
        for(j=0;j<PneigLenght; j++) {
            /* Calculate length edges of face */
            Ec= sqrt(pow2(Pn[j]-P[0])   +pow2(Pn[j+VertexN]-P[1])   +pow2(Pn[j+VertexN*2]-P[2]))+1e-14;
            Eb= sqrt(pow2(Pnop[j]-P[0]) +pow2(Pnop[j+VertexN]-P[1]) +pow2(Pnop[j+VertexN*2]-P[2]))+1e-14;
            Ea= sqrt(pow2(Pn[j]-Pnop[j])+pow2(Pn[j+VertexN]-Pnop[j+VertexN]) +pow2(Pn[j+VertexN*2]-Pnop[j+VertexN*2]))+1e-14;
           
            /* Calculate face surface area */
            s = ((Ea+Eb+Ec)/2);
            h = (2/Ea)*sqrt(s*(s-Ea)*(s-Eb)*(s-Ec))+1e-14;
            x = (pow2(Ea)-pow2(Eb)+pow2(Ec))/(2*Ea);
            
            /* Calculate tangent of 2D triangle */
            Np[0]=-h / sqrt(pow2(-h)+pow2(x));   Np[1]= x / sqrt(pow2(-h)+pow2(x));
            Ns[0]=h / sqrt(pow2(h)+pow2(Ea-x));  Ns[1]=(Ea-x) / sqrt(pow2(h)+pow2(Ea-x));
            Tb[0]= Np[1]+Ns[1] ; Tb[1]=-(Np[0]+Ns[0]);
            
            /* Back to 3D coordinates */
            Pm[0]=(Pn[j]*x+Pnop[j]*(Ea-x))/Ea;
            Pm[1]=(Pn[j+VertexN]*x+Pnop[j+VertexN]*(Ea-x))/Ea;
            Pm[2]=(Pn[j+VertexN*2]*x+Pnop[j+VertexN*2]*(Ea-x))/Ea;
            X3[0]=(Pn[j]-Pnop[j])/Ea;
            X3[1]=(Pn[j+VertexN]-Pnop[j+VertexN])/Ea;
            X3[2]=(Pn[j+VertexN*2]-Pnop[j+VertexN*2])/Ea;
            Y3[0]=(P[0]-Pm[0])/h; 
            Y3[1]=(P[1]-Pm[1])/h;
            Y3[2]=(P[2]-Pm[2])/h;
            
            /* 2D tangent to 3D tangent */
            Tb3D[0]=(X3[0]*Tb[0]+Y3[0]*Tb[1]); 
            Tb3D[1]=(X3[1]*Tb[0]+Y3[1]*Tb[1]); 
            Tb3D[2]=(X3[2]*Tb[0]+Y3[2]*Tb[1]); 
            Tb3D_length=sqrt(pow2(Tb3D[0])+pow2(Tb3D[1])+pow2(Tb3D[2]))+1e-14;
            Tb3D[0]=Tb3D[0]/Tb3D_length;
            Tb3D[1]=Tb3D[1]/Tb3D_length;
            Tb3D[2]=Tb3D[2]/Tb3D_length;
            
            /* Edge Velocity */
            Vv=0.5*(Ec+0.5*Ea);
          
            /* Store the data */
            ETV_index[ETV_num]=i+1;
          
            ETV_index[ETV_num+table_Dims[0]]=Pneig[j];
            
            ET_table[ETV_num]=Tb3D[0];
            ET_table[ETV_num+table_Dims[0]]=Tb3D[1];
            ET_table[ETV_num+2*table_Dims[0]]=Tb3D[2];
            EV_table[ETV_num]=Vv;
            ETV_num++;
            
            
        }

         
    }

    
    /* Remove temporary memory */
    free(Pn);
    free(Pnop);
    
}

