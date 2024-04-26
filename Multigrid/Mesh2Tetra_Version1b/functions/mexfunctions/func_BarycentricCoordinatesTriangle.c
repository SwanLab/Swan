void BarycentricCoordinatesTriangle(double *VN0,double *VN1,double *VN2,double rx, double ry, double *Lambda)
{
    double f12, f20,f01;
    double g12x, g20x, g01x;
    double g12y, g20y, g01y;
    double c12, c20, c01;
    
    /* Normalization factors */
    f12 = ( VN1[1] - VN2[1] ) * VN0[0]  + (VN2[0] - VN1[0] ) * VN0[1] + VN1[0] * VN2[1] - VN2[0] *VN1[1];
    f20 = ( VN2[1] - VN0[1] ) * VN1[0]  + (VN0[0] - VN2[0] ) * VN1[1] + VN2[0] * VN0[1] - VN0[0] *VN2[1];
    f01 = ( VN0[1] - VN1[1] ) * VN2[0]  + (VN1[0] - VN0[0] ) * VN2[1] + VN0[0] * VN1[1] - VN1[0] *VN0[1];

    /* Lambda Gradient */
    g12x = ( VN1[1] - VN2[1] )/f12; g12y = ( VN2[0] - VN1[0] )/f12;
    g20x = ( VN2[1] - VN0[1] )/f20; g20y = ( VN0[0] - VN2[0] )/f20;
    g01x = ( VN0[1] - VN1[1] )/f01; g01y = ( VN1[0] - VN0[0] )/f01;

    /* Center compensation */
    c12 = (VN1[0] * VN2[1] - VN2[0] * VN1[1])/f12;
    c20 = (VN2[0] * VN0[1] - VN0[0] * VN2[1])/f20;
    c01 = (VN0[0] * VN1[1] - VN1[0] * VN0[1])/f01;

    /* Interpolation values */
    Lambda[0]=g12x*rx+g12y*ry+c12;
    Lambda[1]=g20x*rx+g20y*ry+c20;
    Lambda[2]=g01x*rx+g01y*ry+c01;
}
