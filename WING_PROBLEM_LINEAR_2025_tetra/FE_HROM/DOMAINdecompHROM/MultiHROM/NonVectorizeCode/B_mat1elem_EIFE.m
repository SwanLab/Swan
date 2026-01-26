function [Ke,Vall_rot,TRANSF_COORD] = B_mat1elem_EIFE(EIFEoper_all,Xe) ;

if nargin == 0
    load('tmp.mat')
end

ndim = size(Xe,1) ; nnodeE = size(Xe,2)  ;
Ke = zeros(nnodeE*ndim,nnodeE*ndim) ;

% DETERMINATION PARENT DOMAIN 

weig = EIFEoper.INTforces.weights ;  % CECM weights
ngaus = length(weig) ;
nstrain = EIFEoper.MESH.nstrain ;

for  g = 1:ngaus
    istrain = small2large(g,nstrain ) ;
    BmatRED = EIFEoper.INTforces.BmatRED(istrain,:);
    Be = (BmatRED*EIFEoper.OPER.HdefINV_PsiDEFfT*Vall_rot) ;
    celas = EIFEoper.INTforces.MATPRO.celasglo(istrain,:) ;
    Ke = Ke + weig(g)*detJe*(Be'*celas*Be) ;
end