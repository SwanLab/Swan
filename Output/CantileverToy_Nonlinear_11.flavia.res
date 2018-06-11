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

Result "Displacements" "Elastic Problem" 11 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.12919e-03  1.37575e-04 
     3        00000        00000 
     4  1.48899e-03 -4.96750e-04 
End Values

Result "Stress" "Elastic Problem" 11 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.41971e-03  6.38813e-04  4.76250e-05  6.38800e-04 
     2  1.51457e-03  4.75779e-05 -4.76381e-05  4.84138e-04 
End Values

Result "Strain" "Elastic Problem" 11 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1 -5.00000e-01 -5.00000e-01        00000 -5.00000e-01 
     2 -5.00180e-01 -4.99683e-01 -1.79606e-04 -5.00000e-01 
End Values
