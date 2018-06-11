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

Result "Displacements" "Elastic Problem" 5 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  7.00640e-02  1.11602e-02 
     3        00000        00000 
     4  8.77978e-02 -2.72594e-02 
End Values

Result "Stress" "Elastic Problem" 5 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  9.74939e-02  3.84884e-02  4.86315e-03  3.83976e-02 
     2  1.06510e-01  4.59525e-03 -4.97272e-03  2.94805e-02 
End Values

Result "Strain" "Elastic Problem" 5 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  6.32783e-02        00000  5.21531e-03        00000 
     2  7.75035e-02 -4.03959e-02 -5.75266e-03        00000 
End Values
