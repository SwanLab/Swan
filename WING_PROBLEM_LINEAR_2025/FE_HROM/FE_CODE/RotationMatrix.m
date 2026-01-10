function T = RotationMatrix(a)

a = a/180*pi ; 
 
cosTc = (cos(a))^2 ; 
sinTc = (sin(a))^2 ; 
sinTd = (sin(2*a)) ; 
cosTd = cos(2*a) ; 
cosT = (cos(a)) ; 
sinT = (sin(a)) ;

T = [cosTc  sinTc    0      0      0      -sinTd  
  sinTc  cosTc    0      0      0      sinTd   
     0      0       1      0      0        0       
     0      0     0      cosT      sinT           0        
     0      0        0       -sinT      cosT          0         
    0.5*(sinTd)     -0.5*(sinTd)        0       0    0          cosTd         ] ;