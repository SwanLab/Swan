function Gamma = vectTOmatSYMstress(ndim)
if  ndim == 3
    Gamma = sym(zeros(3,3,6));   
    Gamma(1,1,1) = 1 ; Gamma(2,2,2) = 1 ; Gamma(3,3,3) = 1 ;
    Gamma(2,3,4) = 1 ; Gamma(1,3,5) = 1 ; Gamma(1,2,6) = 1 ;
    Gamma(3,2,4) = 1 ; Gamma(3,1,5) = 1 ; Gamma(2,1,6) = 1 ;
else
    Gamma = sym(zeros(2,2,3));  
    Gamma(1,1,1) = 1 ;  Gamma(2,2,2) = 1 ; 
    Gamma(1,2,3) = 1 ;  Gamma(2,1,3) = 1 ;
end