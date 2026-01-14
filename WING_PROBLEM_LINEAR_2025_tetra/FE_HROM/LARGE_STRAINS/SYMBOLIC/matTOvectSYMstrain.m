function GammaBINV = matTOvectSYMstrain(ndim)
if  ndim == 3
    GammaBINV = sym(zeros(6,3,3));   
    GammaBINV(1,1,1) = 1 ;     GammaBINV(2,2,2) = 1 ;     GammaBINV(3,3,3) = 1 ;   
    GammaBINV(4,2,3) = 1 ;     GammaBINV(5,1,3) = 1 ;     GammaBINV(6,1,2) = 1 ;
    GammaBINV(4,3,2) = 1 ;     GammaBINV(5,3,1) = 1 ;     GammaBINV(6,2,1) = 1 ;
     
else
    GammaBINV = sym(zeros(3,2,2));  
    GammaBINV(1,1,1) = 1 ;      GammaBINV(2,2,2) = 1 ;
    GammaBINV(3,1,2) = 1 ;      GammaBINV(3,2,1) = 1 ;
end