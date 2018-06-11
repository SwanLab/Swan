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

Result "Displacements" "Elastic Problem" 15 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.53895e-03  1.87969e-04 
     3        00000        00000 
     4  2.02854e-03 -6.76382e-04 
End Values

Result "Stress" "Elastic Problem" 15 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.93607e-03  8.70450e-04  6.51675e-05  8.70426e-04 
     2  2.06568e-03  6.50796e-05 -6.51920e-05  6.59682e-04 
End Values

Result "Strain" "Elastic Problem" 15 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1 -5.00000e-01 -5.00000e-01        00000 -5.00000e-01 
     2 -5.00244e-01 -4.99568e-01 -2.44294e-04 -5.00000e-01 
End Values
