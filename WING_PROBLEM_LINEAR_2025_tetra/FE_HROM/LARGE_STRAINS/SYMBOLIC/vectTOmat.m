function Lambda = vectTOmat(ndim)
if  ndim == 3
    Lambda = sym(zeros(3,3,9));   
    Lambda(1,1,1) = 1 ;       Lambda(2,2,2) = 1 ;     Lambda(3,3,3) = 1 ;    
    Lambda(2,3,4) = 1;     Lambda(1,3,5) = 1 ;     Lambda(1,2,6) = 1 ;
    Lambda(3,2,7) = 1 ;       Lambda(3,1,8) = 1 ;      Lambda(2,1,9) = 1 ;
else
    Lambda = sym(zeros(2,2,3));  
    Lambda(1,1,1) = 1 ;       Lambda(2,2,2) = 1 ;
    Lambda(1,2,3) = 1 ;          Lambda(2,1,4) = 1 ;
end