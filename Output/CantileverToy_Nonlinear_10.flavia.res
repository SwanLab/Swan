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

Result "Displacements" "Elastic Problem" 10 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.30218e-01  2.29888e-02 
     3        00000        00000 
     4  1.58538e-01 -4.73236e-02 
End Values

Result "Stress" "Elastic Problem" 10 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.95834e-01  6.96307e-02  1.19002e-02  6.92276e-02 
     2  2.19119e-01  1.07871e-02 -1.24703e-02  5.45289e-02 
End Values

Result "Strain" "Elastic Problem" 10 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  1.07856e-01        00000  1.01702e-02        00000 
     2  1.26671e-01 -7.68204e-02 -1.21806e-02        00000 
End Values
