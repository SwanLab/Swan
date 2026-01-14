function Be = QtransfB(BeTILDE,ndim)
%
if nargin == 0
    load('tmp2.mat')
end


nnodeE = size(BeTILDE,2); 
if ndim==2 
   nstrain = 3 ;    Be = zeros(nstrain,nnodeE*ndim) ; 
   column1 = 1:2:(nnodeE*2-1) ;
   column2 = 2:2:nnodeE*2 ;
   Be(1,column1) = BeTILDE(1,:) ; 
   Be(2,column2) = BeTILDE(2,:) ; 
   Be(3,column1) = BeTILDE(2,:) ; 
   Be(3,column2) = BeTILDE(1,:) ; 
elseif ndim ==3 
   nstrain = 6 ;    Be = zeros(nstrain,nnodeE*ndim) ; 
   column1 = 1:3:(nnodeE*3-2) ;
   column2 = 2:3:(nnodeE*3-1) ;
   column3 = 3:3:nnodeE*3 ;    
   Be(1,column1) = BeTILDE(1,:) ; 
   Be(2,column2) = BeTILDE(2,:) ; 
   Be(3,column3) = BeTILDE(3,:) ; 
   Be(4,column2) = BeTILDE(3,:) ; 
   Be(4,column3) = BeTILDE(2,:) ;
   Be(5,column1) = BeTILDE(3,:) ; 
   Be(5,column3) = BeTILDE(1,:) ;
   Be(6,column1) = BeTILDE(2,:) ; 
   Be(6,column2) = BeTILDE(1,:) ;   
else
   error('Incorrect option')
end