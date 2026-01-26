function Be = QtransfBvect_CABLE(BeTILDE,ndim)
%
if nargin == 0
    load('tmp.mat')
end
%Convertion from B-matrix for scalar fielts to B-matrix for vector fields 

nnodeE = size(BeTILDE,2);
nelem = size(BeTILDE,1) ;
 Be = zeros(ndim*nelem,nnodeE*ndim) ;

 for idim = 1:ndim
     indROWS = idim:ndim:(ndim*nelem) ;
     indCOLS = idim:ndim:(ndim*nnodeE) ; 
     
     Be(indROWS,indCOLS) = BeTILDE  ; 
 end


% 
% 
% if ndim==2
%     nstrain = 3 ;    Be = zeros(nstrain*nelem,nnodeE*ndim) ;
%     column1 = 1:2:(nnodeE*2-1) ;
%     column2 = 2:2:nnodeE*2 ;
%     ind1 = 1:nstrain:nstrain*nelem; ind2 = 2:nstrain:nstrain*nelem; ind3 = 3:nstrain:nstrain*nelem ;
%     jnd1 = 1:ndim:ndim*nelem; jnd2 = 2:ndim:ndim*nelem; jnd3 = 3:ndim:ndim*nelem ;
%     Be(ind1,column1) = BeTILDE(jnd1,:) ;
%     Be(ind2,column2) = BeTILDE(jnd2,:) ;
%     Be(ind3,column1) = BeTILDE(jnd2,:) ;
%     Be(ind3,column2) = BeTILDE(jnd1,:) ;
% elseif ndim ==3
%     nstrain = 6 ;    Be = zeros(nstrain*nelem,nnodeE*ndim) ;
%     column1 = 1:3:(nnodeE*3-2) ;
%     column2 = 2:3:(nnodeE*3-1) ;
%     column3 = 3:3:nnodeE*3 ;
%     ind1 = 1:nstrain:nstrain*nelem; ind2 = 2:nstrain:nstrain*nelem; ind3 = 3:nstrain:nstrain*nelem ;
%     ind4 = 4:nstrain:nstrain*nelem; ind5 = 5:nstrain:nstrain*nelem; ind6 = 6:nstrain:nstrain*nelem ;
%     jnd1 = 1:ndim:ndim*nelem; jnd2 = 2:ndim:ndim*nelem; jnd3 = 3:ndim:ndim*nelem ;
%     jnd4 = 4:ndim:ndim*nelem; jnd5 = 5:ndim:ndim*nelem; jnd6 = 6:ndim:ndim*nelem ;
%     Be(ind1,column1) = BeTILDE(jnd1,:) ;
%     Be(ind2,column2) = BeTILDE(jnd2,:) ;
%     Be(ind3,column3) = BeTILDE(jnd3,:) ;
%     Be(ind4,column2) = BeTILDE(jnd3,:) ;
%     Be(ind4,column3) = BeTILDE(jnd2,:) ;
%     Be(ind5,column1) = BeTILDE(jnd3,:) ;
%     Be(ind5,column3) = BeTILDE(jnd1,:) ;
%     Be(ind6,column1) = BeTILDE(jnd2,:) ;
%     Be(ind6,column2) = BeTILDE(jnd1,:) ;
% else
%     error('Incorrect option')
% end