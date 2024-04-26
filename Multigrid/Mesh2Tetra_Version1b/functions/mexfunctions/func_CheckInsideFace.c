
bool CheckInsideFace(double *VN0, double *VN1, double *VN2, double rx, double ry, double *mi, double *ma, bool ignore_corners) {
    bool inter;
    double Lambda[3];
    double mv, vv;
    if((rx<mi[0])||(rx>ma[0])||(ry<mi[1])||(ry>ma[1]))
    {
        inter=false;
        return inter;
    }

    BarycentricCoordinatesTriangle(VN0, VN1, VN2, rx, ry, Lambda);
    
    /* Check if inside face */
    if(ignore_corners) {
        mv=1e-8; vv=1-1e-8;
    }
    else {
        mv=0; vv=1;
    }
    
    inter=(Lambda[0]>=mv)&&(Lambda[0]<=vv)&&(Lambda[1]>=mv)&&(Lambda[1]<=vv)&& (Lambda[2]>=mv)&&(Lambda[2]<=vv);
    return inter;
}
