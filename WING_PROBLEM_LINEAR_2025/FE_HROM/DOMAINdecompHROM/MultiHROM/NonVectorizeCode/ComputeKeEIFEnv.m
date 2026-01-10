function [Ke,Vall_rot,TRANSF_COORD] = ComputeKeEIFEnv(EIFEoper,BeALL,weig) ;

if nargin == 0
    load('tmp.mat')
end

%ndim = size(Xe,1) ; nnodeE = size(Xe,2)  ;
Ke = zeros(size(BeALL,2)) ;


ngaus = length(weig) ;
nstrain = EIFEoper.MESH.nstrain ;

for  g = 1:ngaus
    istrain = small2large(g,nstrain ) ;
   % BmatRED = EIFEoper.INTforces.BmatRED(istrain,:);
    Be = BeALL(istrain,:) ;
    celas = EIFEoper.INTforces.MATPRO.celasglo(istrain,:) ;
    Ke = Ke + weig(g)*(Be'*celas*Be) ;
end