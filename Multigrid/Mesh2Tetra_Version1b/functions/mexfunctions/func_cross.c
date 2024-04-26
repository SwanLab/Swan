void cross(double *a, double *b, double *n) {
    n[0]=a[1]*b[2]-a[2]*b[1]; n[1]=a[2]*b[0]-a[0]*b[2]; n[2]=a[0]*b[1]-a[1]*b[0];
}
