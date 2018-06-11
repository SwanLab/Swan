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

Result "Displacements" "Elastic Problem" 3 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  4.35775e-02  6.43764e-03 
     3        00000        00000 
     4  5.55414e-02 -1.76551e-02 
End Values

Result "Stress" "Elastic Problem" 3 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  5.83575e-02  2.41949e-02  2.57875e-03  2.41653e-02 
     2  6.31322e-02  2.48675e-03 -2.61193e-03  1.84104e-02 
End Values

Result "Strain" "Elastic Problem" 3 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  4.08668e-02        00000  3.08482e-03        00000 
     2  5.12702e-02 -2.48437e-02 -3.27814e-03        00000 
End Values
