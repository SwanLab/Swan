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

Result "Displacements" "Elastic Problem" 1 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.51034e-02  1.98874e-03 
     3        00000        00000 
     4  1.96701e-02 -6.44645e-03 
End Values

Result "Stress" "Elastic Problem" 1 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.93935e-02  8.49028e-03  7.23314e-04  8.48752e-03 
     2  2.07785e-02  7.13951e-04 -7.26179e-04  6.43597e-03 
End Values

Result "Strain" "Elastic Problem" 1 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  1.47661e-02        00000  9.79574e-04        00000 
     2  1.91123e-02 -8.52371e-03 -1.00021e-03        00000 
End Values
