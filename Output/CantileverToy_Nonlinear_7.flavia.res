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

Result "Displacements" "Elastic Problem" 7 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  9.50790e-02  1.59321e-02 
     3        00000        00000 
     4  1.17589e-01 -3.58565e-02 
End Values

Result "Stress" "Elastic Problem" 7 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.36757e-01  5.16366e-02  7.48202e-03  5.14477e-02 
     2  1.50863e-01  6.94385e-03 -7.72576e-03  3.98775e-02 
End Values

Result "Strain" "Elastic Problem" 7 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  8.27939e-02        00000  7.27750e-03        00000 
     2  9.94764e-02 -5.53380e-02 -8.31062e-03        00000 
End Values
