void SphereFrom4Points(double *PA, double *PB, double *PC, double *PD, double *Cx, double *Cy, double *Cz, double *L) {
    double x1, x2, x3, x4;
    double **a;
    double m11, m12, m13, m14, m15;
    int i;
    
    a = malloc(4*sizeof(double *)); for (i=0;i<4;i++) { a[i] = malloc(4*sizeof(double)); }
    
    x1=PA[0]*PA[0]+PA[1]*PA[1]+PA[2]*PA[2];
    x2=PB[0]*PB[0]+PB[1]*PB[1]+PB[2]*PB[2];
    x3=PC[0]*PC[0]+PC[1]*PC[1]+PC[2]*PC[2];
    x4=PD[0]*PD[0]+PD[1]*PD[1]+PD[2]*PD[2];
    
    a[0][0]=PA[0]; a[0][1]=PA[1]; a[0][2]=PA[2]; a[0][3]=1;
    a[1][0]=PB[0]; a[1][1]=PB[1]; a[1][2]=PB[2]; a[1][3]=1;
    a[2][0]=PC[0]; a[2][1]=PC[1]; a[2][2]=PC[2]; a[2][3]=1;
    a[3][0]=PD[0]; a[3][1]=PD[1]; a[3][2]=PD[2]; a[3][3]=1;
    m11 = det(a);
    
    a[0][0]=x1; a[0][1]=PA[1]; a[0][2]=PA[2]; a[0][3]=1;
    a[1][0]=x2; a[1][1]=PB[1]; a[1][2]=PB[2]; a[1][3]=1;
    a[2][0]=x3; a[2][1]=PC[1]; a[2][2]=PC[2]; a[2][3]=1;
    a[3][0]=x4; a[3][1]=PD[1]; a[3][2]=PD[2]; a[3][3]=1;
    m12 = det(a);
    
    a[0][0]=PA[0]; a[0][1]= x1; a[0][2]=PA[2]; a[0][3]=1;
    a[1][0]=PB[0]; a[1][1]= x2; a[1][2]=PB[2]; a[1][3]=1;
    a[2][0]=PC[0]; a[2][1]= x3; a[2][2]=PC[2]; a[2][3]=1;
    a[3][0]=PD[0]; a[3][1]= x4; a[3][2]=PD[2]; a[3][3]=1;
    m13 = det(a);
    
    a[0][0]=PA[0]; a[0][1]=PA[1]; a[0][2]=x1; a[0][3]=1;
    a[1][0]=PB[0]; a[1][1]=PB[1]; a[1][2]=x2; a[1][3]=1;
    a[2][0]=PC[0]; a[2][1]=PC[1]; a[2][2]=x3; a[2][3]=1;
    a[3][0]=PD[0]; a[3][1]=PD[1]; a[3][2]=x4; a[3][3]=1;
    m14 = det(a);
    
    a[0][0]=x1; a[0][1]=PA[0]; a[0][2]=PA[1]; a[0][3]=PA[2];
    a[1][0]=x2; a[1][1]=PB[0]; a[1][2]=PB[1]; a[1][3]=PB[2];
    a[2][0]=x3; a[2][1]=PC[0]; a[2][2]=PC[1]; a[2][3]=PC[2];
    a[3][0]=x4; a[3][1]=PD[0]; a[3][2]=PD[1]; a[3][3]=PD[2];
    m15 = det(a);
    
    if (m11 == 0) { printf("Points define no sphere \n"); }
    
    Cx[0] = 0.5 * m12 / m11;
    Cy[0] = 0.5 * m13 / m11;
    Cz[0] = 0.5 * m14 / m11;
    L[0] = sqrt(Cx[0]*Cx[0] + Cy[0]*Cy[0] + Cz[0]*Cz[0] - m15/m11);
    
    for (i=0;i<4;i++) { free(a[i]); } free(a);
    
}
