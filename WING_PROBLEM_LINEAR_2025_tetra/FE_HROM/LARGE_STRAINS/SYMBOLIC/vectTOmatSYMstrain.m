function GammaB = vectTOmatSYMstrain(ndim)
if  ndim == 3
    GammaB = sym(zeros(3,3,6));   
    GammaB(1,1,1) = 1 ;     GammaB(2,2,2) = 1 ;     GammaB(3,3,3) = 1 ;  
    GammaB(2,3,4) = 0.5 ;   GammaB(1,3,5) = 0.5 ;   GammaB(1,2,6) = 0.5 ;
    GammaB(3,2,4) = 0.5 ;   GammaB(3,1,5) = 0.5 ;   GammaB(2,1,6) = 0.5 ;
else
    GammaB = sym(zeros(2,2,3));  
    GammaB(1,1,1) = 1 ;      GammaB(2,2,2) = 1 ;
    GammaB(1,2,3) = 0.5 ;        GammaB(2,1,3) = 0.5 ;
end