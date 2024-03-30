GiD Post Results File 1.0 
GaussPoints "Gauss Points" Elemtype Quadrilateral
Number of Gauss Points: 4 
Nodes not included 
Natural Coordinates: given 
-5.77350e-01 -5.77350e-01 
 5.77350e-01 -5.77350e-01 
-5.77350e-01  5.77350e-01 
 5.77350e-01  5.77350e-01 
End GaussPoints  

Result "displacement" "FunResults" 0 Vector OnNodes 
ComponentNames  "x1", "x2"
Values  
1        00000        00000 
2 -3.11111e-02        00000 
3 -3.11111e-02  9.33333e-02 
4  2.07801e-17  9.33333e-02 
End Values


Result "strain" "FunResults" 0 Vector OnGaussPoints "Gauss Points" 
ComponentNames  "x1", "x2", "x3"
Values  
1 -1.55556e-02  4.66667e-02  2.59293e-18 
 -1.55556e-02  4.66667e-02  6.06238e-18 
 -1.55556e-02  4.66667e-02  1.17097e-18 
 -1.55556e-02  4.66667e-02  9.43590e-18 
 -1.55556e-02  4.66667e-02  6.83381e-18 
 -1.55556e-02  4.66667e-02  4.64042e-18 
 -1.55556e-02  4.66667e-02  7.24251e-18 
 -1.55556e-02  4.66667e-02  3.36436e-18 
 -1.55556e-02  4.66667e-02  8.66446e-18 
End Values

