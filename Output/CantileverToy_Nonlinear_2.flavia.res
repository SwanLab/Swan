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

Result "Displacements" "Elastic Problem" 2 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  2.96116e-02  4.15762e-03 
     3        00000        00000 
     4  3.81239e-02 -1.22918e-02 
End Values

Result "Stress" "Elastic Problem" 2 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  3.88500e-02  1.65390e-02  1.58980e-03  1.65268e-02 
     2  4.18244e-02  1.55040e-03 -1.60297e-03  1.25553e-02 
End Values

Result "Strain" "Elastic Problem" 2 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  2.83383e-02        00000  2.01929e-03        00000 
     2  3.60713e-02 -1.67928e-02 -2.10408e-03        00000 
End Values
