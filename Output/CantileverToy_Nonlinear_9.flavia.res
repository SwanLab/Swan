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

Result "Displacements" "Elastic Problem" 9 Vector OnNodes ""
ComponentNames  "Ux", "Uy"
Values
     1        00000        00000 
     2  1.18834e-01  2.06582e-02 
     3        00000        00000 
     4  1.45394e-01 -4.36919e-02 
End Values

Result "Stress" "Elastic Problem" 9 Vector OnGaussPoints "Guass up?"
ComponentNames  "Sx", "Sy", "Sxy", "Sz"
Values
     1  1.76120e-01  6.38372e-02  1.03704e-02  6.35142e-02 
     2  1.96144e-01  9.47071e-03 -1.08141e-02  4.97528e-02 
End Values

Result "Strain" "Elastic Problem" 9 Vector OnGaussPoints "Guass up?"
ComponentNames  "Ex", "Ey", "Exy", "Ez"
Values
     1  9.99900e-02        00000  9.23498e-03        00000 
     2  1.18250e-01 -6.97716e-02 -1.08922e-02        00000 
End Values
