/*
 * Recursive definition of determinate using expansion by minors.
 */
double Determinant(double **a, int n) {
    int i, j, j1, j2;
    double det = 0;
    double **m;
    if (n < 1) {
    } else if (n == 1) { det = a[0][0];
    } else if (n == 2) { det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else { det = 0;
    for (j1=0;j1<n;j1++) {
        m = malloc((n-1)*sizeof(double *));
        for (i=0;i<n-1;i++) m[i] = malloc((n-1)*sizeof(double));
        for (i=1;i<n;i++) {
            j2 = 0; for (j=0;j<n;j++) { if (j == j1) {continue; } m[i-1][j2] = a[i][j]; j2++; }
        }
        det += pow(-1.0, 1.0+j1+1.0) * a[0][j1] * Determinant(m, n-1);
        for (i=0;i<n-1;i++) { free(m[i]); }  free(m);
    }
    }
    return(det);
}

__inline double det(double **a) { return Determinant(a, 4); }
