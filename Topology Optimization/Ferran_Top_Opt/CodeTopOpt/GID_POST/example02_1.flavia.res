Gid Post Results File 1.0 
### 
# ENG_AERO_COMP  V.1.0 
# 
GaussPoints "My Gauss" ElemType Number Of Gauss Points: 1
Natural Coordinates: Internal
End gausspoints
Result "vel" "time"       00001 Vector OnNodes 
Values 
     1         00000         00000 
     2         00000         00000 
     3         00000 -9.09923e-001 
     4         00000 -5.52628e-001 
     5         00000        -00001 
     6         00000        -00002 
End Values 
# 
Result "pressure" "time"       00001 Scalar  OnGaussPoints "My Gauss"
Values 
     1 -2.00475e-001 
     2 -3.25376e-001 
     3 -3.25704e-001 
     4 -2.42826e-001 
     5 6.46427e-003 
     6 -3.40210e-001 
End Values 
# 
