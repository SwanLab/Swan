function [BBnw,BB,weigEM,weigEMdiag,posgp] = BMatQuad_JAHO(COOR,CONNECT,BBAR,DATA)

if nargin == 0
    load('tmp0.mat')
end

 posgp = 1/sqrt(3)*[-1 1 1 -1
            -1 -1 1 1 ]; 

[N dN_dXi] = QuadElem_N_derN  ;
weightsE = [1 1 1 1] ;
nelem = size(CONNECT,1) ;
ngausE = size(N,1) ;
ndime = size(COOR,2) ;
nnode = size(COOR,1) ;
nstrain  =4 ;
nnodeE = size(N,2) ;
ngaus = nelem*nstrain ;

BBnwNA = zeros(nstrain*ngausE*nelem,nnodeE*ndime);
weigEM = zeros(ngausE*nelem,1) ;
disp(['Constructing element B-matrices...'])
for ielem = 1:nelem
    nodes = CONNECT(ielem,:) ;
    X = COOR(nodes,:);
    if BBAR==0
        [BBnwLOC detJ] = BmatQuad_loc(X,dN_dXi);
    else
        [BBnwLOC detJ] = BmatQuad_locBBAR(X,dN_dXi,weightsE);
    end
    iniROW = (ielem-1)*nstrain*ngausE +1; finROW = ielem*nstrain*ngausE ;
    BBnwNA(iniROW:finROW,:) = BBnwLOC ;
    weigL = detJ.*weightsE ;
    iniROW = (ielem-1)*ngausE +1; finROW = ielem*ngausE ;
    weigEM(iniROW:finROW) = weigL' ;
end

%%% Assembling BBnw
% -----------------
disp(['Assembling   BBnw...']) ;

BBnw = AsseblyBBnw(BBnwNA,nstrain,nelem,nnodeE,ndime,ngausE,CONNECT,nnode) ;

weigEMaux =   repmat(weigEM',nstrain,1) ;
weigEMdiag = reshape(weigEMaux,nstrain*ngaus,1);
weigEMdiag = sparse(1:ngaus*nstrain,1:ngaus*nstrain,weigEMdiag);
disp(['Calculating   BB...']) ;


if DATA.ComputeBstW == 1
BB = weigEMdiag*BBnw ;
else
   BB = [] ;  
   weigEMdiag = diag(weigEMdiag) ; 
end
%%%
% TESTIMP = 1 ;
% if TESTIMP == 1
%     % TEnsile deformation
%     u = 0.2 ;
%     nu = 0.3 ;
%     d = [0 0 u 0  u  u*nu  0 u*nu]' ;
%     BBnw*d;
%     % Shear deformation
%     d = [0 0 0 0  u  0 u  0]' ;
%     BBnw*d
% end