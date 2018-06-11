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

Result "Displacements" "Elastic Problem" 6 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  8.26626e-02  1.35494e-02 
     3        00000        00000 
     4  1.02838e-01 -3.16134e-02 
End Values

Result "Stress" "Elastic Problem" 6 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.17112e-01  4.51920e-02  6.13533e-03  4.50567e-02 
     2  1.28568e-01  5.74398e-03 -6.30421e-03  3.47524e-02 
End Values

Result "Strain" "Elastic Problem" 6 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  7.33580e-02        00000  6.25803e-03        00000 
     2  8.89460e-02 -4.79360e-02 -7.02578e-03        00000 
End Values
