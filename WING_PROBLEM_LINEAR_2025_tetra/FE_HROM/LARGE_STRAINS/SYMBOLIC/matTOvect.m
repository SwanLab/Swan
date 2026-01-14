function LambdaINV = matTOvect(ndim)
if  ndim == 3
    LambdaINV = sym(zeros(9,3,3));   
    LambdaINV(1,1,1) = 1 ;      LambdaINV(2,2,2) = 1 ;
    LambdaINV(3,3,3) = 1 ;       LambdaINV(4,2,3) = 1 ;
    LambdaINV(5,1,3) = 1 ;       LambdaINV(6,1,2) = 1 ;
    LambdaINV(7,3,2) = 1 ;       LambdaINV(8,3,1) = 1 ;
    LambdaINV(9,2,1) = 1 ;
     
else
    LambdaINV = sym(zeros(4,2,2));  
    LambdaINV(1,1,1) = 1 ;      LambdaINV(2,2,2) = 1 ;
    LambdaINV(3,1,2) = 1 ;        LambdaINV(4,2,1) = 1 ;
end