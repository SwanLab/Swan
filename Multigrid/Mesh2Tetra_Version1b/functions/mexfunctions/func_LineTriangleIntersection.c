bool LineTriangleIntersection(double *A, double *B, double *C, double *P1, double *P2, bool ignore_corners, double *inter_xyz) {
    bool intersect;
    double T1[3], T2[3], N[3];
    double l;
    double Pline[3], Vline[3], P[3];
    double P_2D[2], A_2D[2], B_2D[2], C_2D[2];
    double mi[2], ma[2];
    int i, j;
    double t;
    bool check;    
    check=(mind(mind(A[0],B[0]),C[0])>maxd(P1[0],P2[0]))||(mind(mind(A[1],B[1]),C[1])>maxd(P1[1],P2[1]))||(mind(mind(A[2],B[2]),C[2])>maxd(P1[2],P2[2]))||(maxd(maxd(A[0],B[0]),C[0])<mind(P1[0],P2[0]))||(maxd(maxd(A[1],B[1]),C[1])<mind(P1[1],P2[1]))||(maxd(maxd(A[2],B[2]),C[2])<mind(P1[2],P2[2]));
    if(check) { return false; }

    
    /* http:/*www.cs.brown.edu/~scd/facts.html */
    
    /* Calculate normal of triangle */
    for(i=0; i<3; i++){
        T1[i]=A[i]-C[i]; T2[i]=B[i]-C[i];
    }
    /* N = cross(A-C,B-C); */
    N[0]=T1[1]*T2[2]-T1[2]*T2[1];
    N[1]=T1[2]*T2[0]-T1[0]*T2[2];
    N[2]=T1[0]*T2[1]-T1[1]*T2[0];
    /* N=N./sqrt(sum(N.^2)); */
    l=sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    N[0]/=l; N[1]/=l; N[2]/=l;
    
    for(i=0; i<3; i++){
        Pline[i] = P1[i]; Vline[i] = P2[i]-P1[i]; T1[i]=C[i]-Pline[i];
    }
    
    /* Normalized plane intersection position on line [0..1] */
    t = (N[0]*T1[0]+N[1]*T1[1]+N[2]*T1[2]) / ((N[0]*Vline[0]+N[1]*Vline[1]+N[2]*Vline[2])+1e-16);
    
    /* If no intersection between points defining the line, return */
    if((t<0)||(t>1)) {
        intersect=false;
        return intersect;
    }
    
    /* 3D xyz Intersection Point */
    for(i=0; i<3; i++){
        P[i] =  Pline[i] + t * Vline[i];
    }
	inter_xyz[0]=P[0];
	inter_xyz[1]=P[1];
	inter_xyz[2]=P[2];
	
    
    /* Drop x,y or z to get a 2D triangle intersection problem */
    if(N[0]>N[1]) {
        if(N[0]>N[2]) {
            i=1; j=2;
        }
        else {
            i=0; j=1;
        }
    }
    else {
        if(N[1]>N[2]) {
            i=0; j=2;
        }
        else {
            i=0; j=1;
        }
    }
    
    P_2D[0]=P[i]; P_2D[1]=P[j];
    A_2D[0]=A[i]; A_2D[1]=A[j];
    B_2D[0]=B[i]; B_2D[1]=B[j];
    C_2D[0]=C[i]; C_2D[1]=C[j];
    
    /* Boundary box */
    mi[0]=mind(mind(A_2D[0],B_2D[0]),C_2D[0]);
    mi[1]=mind(mind(A_2D[1],B_2D[1]),C_2D[1]);
    
    ma[0]=maxd(maxd(A_2D[0],B_2D[0]),C_2D[0]);
    ma[1]=maxd(maxd(A_2D[1],B_2D[1]),C_2D[1]);
    
    /* Check if the insection point is inside the face */
    intersect=CheckInsideFace(A_2D, B_2D, C_2D, P_2D[0], P_2D[1], mi, ma,  ignore_corners);
    return intersect;
}
