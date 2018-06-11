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

Result "Displacements" "Elastic Problem" 12 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.23168e-03  1.50156e-04 
     3        00000        00000 
     4  1.62397e-03 -5.41708e-04 
End Values

Result "Stress" "Elastic Problem" 12 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.54880e-03  6.96755e-04  5.19995e-05  6.96740e-04 
     2  1.65233e-03  5.19433e-05 -5.20151e-05  5.28050e-04 
End Values

Result "Strain" "Elastic Problem" 12 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1 -5.00000e-01 -5.00000e-01        00000 -5.00000e-01 
     2 -5.00196e-01 -4.99654e-01 -1.95809e-04 -5.00000e-01 
End Values
