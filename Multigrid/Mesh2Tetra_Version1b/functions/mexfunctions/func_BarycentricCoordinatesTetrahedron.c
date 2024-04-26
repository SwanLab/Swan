void BarycentricCoordinatesTetrahedron(double *p1, double *p2, double *p3, double *p4, double rx, double ry, double rz, double *Lambda) {
    int i;
    double rt[3];
    double r1[3], r2[3], r3[3];
    double r4[3], r5[3], r6[3];
    double g1[3], g2[3], g3[3], g4[3];
    double c1, c2, c3, c4;
    double r_b[3];
    double r[3];
    double J;
    
    /* Edge Vectors */
    for(i=0; i<3; i++) {
        r1[i] = p2[i]-p1[i];  r2[i] = p3[i]-p1[i];
        r3[i] = p4[i]-p1[i];  r4[i] = p3[i]-p2[i];
        r5[i] = p2[i]-p4[i];  r6[i] = p4[i]-p3[i];
    }
    
    /* Jacobian determinant */
    cross(r1, r2, rt);
    J = rt[0]*r3[0]+rt[1]*r3[1]+rt[2]*r3[2];
    
    /* Lambda Gradient */
    cross(r4, r5, rt); g1[0] = rt[0]/J; g1[1] = rt[1]/J; g1[2] = rt[2]/J;
    cross(r2, r3, rt); g2[0] = rt[0]/J; g2[1] = rt[1]/J; g2[2] = rt[2]/J;
    cross(r3, r1, rt); g3[0] = rt[0]/J; g3[1] = rt[1]/J; g3[2] = rt[2]/J;
    cross(r1, r2, rt); g4[0] = rt[0]/J; g4[1] = rt[1]/J; g4[2] = rt[2]/J;
    
    /* Center compensation */
    for(i=0; i<3; i++) {
        r_b[i] = -(p1[i] + p2[i] + p3[i] + p4[i])/4.0;
    }
    c1=(r_b[0]*g1[0]+r_b[1]*g1[1]+r_b[2]*g1[2])+1/4.0;
    c2=(r_b[0]*g2[0]+r_b[1]*g2[1]+r_b[2]*g2[2])+1/4.0;
    c3=(r_b[0]*g3[0]+r_b[1]*g3[1]+r_b[2]*g3[2])+1/4.0;
    c4=(r_b[0]*g4[0]+r_b[1]*g4[1]+r_b[2]*g4[2])+1/4.0;
    
    /* Current location */
    r[0]=rx;
    r[1]=ry;
    r[2]=rz;
    
    /* Interpolation values */
    Lambda[0]=r[0]*g1[0]+r[1]*g1[1]+r[2]*g1[2]+c1;
    Lambda[1]=r[0]*g2[0]+r[1]*g2[1]+r[2]*g2[2]+c2;
    Lambda[2]=r[0]*g3[0]+r[1]*g3[1]+r[2]*g3[2]+c3;
    Lambda[3]=r[0]*g4[0]+r[1]*g4[1]+r[2]*g4[2]+c4;
}
