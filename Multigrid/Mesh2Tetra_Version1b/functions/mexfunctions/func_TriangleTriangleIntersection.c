bool TriangleTriangleIntersection(double *P1,double *P2,double *P3,double *O1,double *O2,double *O3, bool ignore_corners)
{
    int i=0;
    double A[3],B[3],C[3], E1[3], E2[3];
	double inter_xyz[3];
    bool check;    
    check=(mind(mind(P1[0],P2[0]),P3[0])>maxd(maxd(O1[0],O2[0]),O3[0]))||(mind(mind(P1[1],P2[1]),P3[1])>maxd(maxd(O1[1],O2[1]),O3[1]))||(mind(mind(P1[2],P2[2]),P3[2])>maxd(maxd(O1[2],O2[2]),O3[2]))||(maxd(maxd(P1[0],P2[0]),P3[0])<mind(mind(O1[0],O2[0]),O3[0]))||(maxd(maxd(P1[1],P2[1]),P3[1])<mind(mind(O1[1],O2[1]),O3[1]))||(maxd(maxd(P1[2],P2[2]),P3[2])<mind(mind(O1[2],O2[2]),O3[2]));
    if(check) { return false; }
    
   
    for(i=0; i<3; i++){
        A[i]=P1[i]; B[i]=P2[i]; C[i]=P3[i]; E1[i]=O1[i]; E2[i]=O2[i];
    }
    if(LineTriangleIntersection(A,B,C,E1,E2,ignore_corners,inter_xyz)) { return true; }

    for(i=0; i<3; i++){
        A[i]=P1[i]; B[i]=P2[i]; C[i]=P3[i]; E1[i]=O2[i]; E2[i]=O3[i];
    }
    if(LineTriangleIntersection(A,B,C,E1,E2,ignore_corners,inter_xyz)) { return true; }

    for(i=0; i<3; i++){
        A[i]=P1[i]; B[i]=P2[i]; C[i]=P3[i]; E1[i]=O3[i]; E2[i]=O1[i];
    }
    if(LineTriangleIntersection(A,B,C,E1,E2,ignore_corners,inter_xyz)) { return true; }

    for(i=0; i<3; i++){
        A[i]=O1[i]; B[i]=O2[i]; C[i]=O3[i]; E1[i]=P1[i]; E2[i]=P2[i];
    }
    if(LineTriangleIntersection(A,B,C,E1,E2,ignore_corners,inter_xyz)) { return true; }

    for(i=0; i<3; i++){
        A[i]=O1[i]; B[i]=O2[i]; C[i]=O3[i]; E1[i]=P2[i]; E2[i]=P3[i];
    }
    if(LineTriangleIntersection(A,B,C,E1,E2,ignore_corners,inter_xyz)) { return true; }

    for(i=0; i<3; i++){
        A[i]=O1[i]; B[i]=O2[i]; C[i]=O3[i]; E1[i]=P3[i]; E2[i]=P1[i];
    }
    if(LineTriangleIntersection(A,B,C,E1,E2,ignore_corners,inter_xyz)) { return true; }
    return false;
}
