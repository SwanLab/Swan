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

Result "Displacements" "Elastic Problem" 14 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.43655e-03  1.75352e-04 
     3        00000        00000 
     4  1.89374e-03 -6.31524e-04 
End Values

Result "Stress" "Elastic Problem" 14 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.80698e-03  8.12574e-04  6.07707e-05  8.12553e-04 
     2  1.92788e-03  6.06942e-05 -6.07921e-05  6.15821e-04 
End Values

Result "Strain" "Elastic Problem" 14 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1 -5.00000e-01 -5.00000e-01        00000 -5.00000e-01 
     2 -5.00228e-01 -4.99597e-01 -2.28153e-04 -5.00000e-01 
End Values
