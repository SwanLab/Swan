clc
clear all

syms     F1 F2 F3 F4 F5 F6 F7 F8 F9  S1 S2 S3 S4 S5 S6
ndim =2 ;
if ndim == 3
    F = [F1  F6 F5
        F9  F2 F4
        F8   F7 F3 ] ;
    
    S  = [S1 S6 S5
        S6 S2 S4
        S5 S4 S3] ;
else
    F = [F1  F3
        F4  F2  ] ;
    
    S  = [S1 S3; 
          S3 S2] ;
    
end

P = F*S

if ndim == 2
      P1 = P(1,1)
    P2 = P(2,2)
     P3 = P(1,2)
    P4 = P(2,1)
else
    P1 = P(1,1)
    P2 = P(2,2)
    P3 = P(3,3)
    P4 = P(2,3)
    P5 = P(1,3)
    P6 = P(1,2)
    P7 = P(3,2)
    P8 = P(3,1)
    P9 = P(2,1)
end