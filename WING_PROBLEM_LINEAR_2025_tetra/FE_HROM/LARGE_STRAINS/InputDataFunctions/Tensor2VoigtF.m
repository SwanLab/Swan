function Fvoigt = Tensor2VoigtF(Ftensor) 

ndim = size(Ftensor,1) ; 
ndimNEW = ndim^2; 

Fvoigt = zeros(ndimNEW,1) ; 


if ndim == 3
    
   Fvoigt(1) = Ftensor(1,1) ; 
   Fvoigt(2) = Ftensor(2,2) ;
   Fvoigt(3) = Ftensor(3,3) ;
   Fvoigt(4) = Ftensor(2,3) ;
   Fvoigt(5) = Ftensor(1,3) ;
   Fvoigt(6) = Ftensor(1,2) ;
   Fvoigt(7) = Ftensor(3,2) ;
   Fvoigt(8) = Ftensor(3,1) ;
   Fvoigt(9) = Ftensor(2,1) ;
else
   Fvoigt(1) = Ftensor(1,1) ; 
   Fvoigt(2) = Ftensor(2,2) ;
   Fvoigt(3) = Ftensor(1,2) ;
   Fvoigt(4) = Ftensor(2,1) ;
    
end