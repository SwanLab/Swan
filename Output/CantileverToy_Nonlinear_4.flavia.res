GiD Post Results File 1.0

####################################################
################# FEM-MAT-OO v.1.0 #################
####################################################

GaussPoints "Guass up?" Elemtype Triangle
Number of Gauss Points: 1
Nodes not included
Natural Coordinates: given
 3.33333e-01  3.33333e-01 
End GaussPoints

Result "Displacements" "Elastic Problem" 4 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  5.70478e-02  8.78236e-03 
     3        00000        00000 
     4  7.20623e-02 -2.26229e-02 
End Values

Result "Stress" "Elastic Problem" 4 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  7.79076e-02  3.14996e-02  3.67413e-03  3.14439e-02 
     2  8.46961e-02  3.50606e-03 -3.73892e-03  2.40417e-02 
End Values

Result "Strain" "Elastic Problem" 4 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  5.24780e-02        00000  4.15470e-03        00000 
     2  6.50049e-02 -3.27043e-02 -4.49964e-03        00000 
End Values
