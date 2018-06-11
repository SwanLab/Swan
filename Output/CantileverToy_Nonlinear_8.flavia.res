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

Result "Displacements" "Elastic Problem" 8 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.07129e-01  1.83037e-02 
     3        00000        00000 
     4  1.31761e-01 -3.98769e-02 
End Values

Result "Stress" "Elastic Problem" 8 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.56428e-01  5.78450e-02  8.89582e-03  5.75935e-02 
     2  1.73390e-01  8.18800e-03 -9.23090e-03  4.48730e-02 
End Values

Result "Strain" "Elastic Problem" 8 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  9.16523e-02        00000  8.27044e-03        00000 
     2  1.09212e-01 -6.26134e-02 -9.60098e-03        00000 
End Values
