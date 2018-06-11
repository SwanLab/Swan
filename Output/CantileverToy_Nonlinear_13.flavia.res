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

Result "Displacements" "Elastic Problem" 13 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.33413e-03  1.62748e-04 
     3        00000        00000 
     4  1.75889e-03 -5.86633e-04 
End Values

Result "Stress" "Elastic Problem" 13 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.67788e-03  7.54675e-04  5.63814e-05  7.54657e-04 
     2  1.79010e-03  5.63155e-05 -5.63998e-05  5.71944e-04 
End Values

Result "Strain" "Elastic Problem" 13 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1 -5.00000e-01 -5.00000e-01        00000 -5.00000e-01 
     2 -5.00212e-01 -4.99625e-01 -2.11991e-04 -5.00000e-01 
End Values
