function GammaINV = matTOvectSYMstress(ndim)
if  ndim == 3
    GammaINV = sym(zeros(6,3,3));   
    GammaINV(1,1,1) = 1 ;     GammaINV(2,2,2) = 1 ;     GammaINV(3,3,3) = 1 ;
    GammaINV(4,2,3) = 0.5 ;   GammaINV(5,1,3) = 0.5 ;   GammaINV(6,1,2) = 0.5 ;    
    GammaINV(4,3,2) = 0.5 ;   GammaINV(5,3,1) = 0.5 ;   GammaINV(6,2,1) = 0.5 ;
else
    GammaINV = sym(zeros(3,2,2));  
    GammaINV(1,1,1) = 1 ;     GammaINV(2,2,2) = 1 ;
    GammaINV(3,1,2) = 0.5 ;    GammaINV(3,2,1) = 0.5 ;    
end