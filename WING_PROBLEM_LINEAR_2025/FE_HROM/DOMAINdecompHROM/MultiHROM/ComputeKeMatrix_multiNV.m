function Ke = ComputeKeMatrix_multiNV(EIFEoper,Xe) ;

if nargin == 0
    load('tmp1.mat')
end

weig = EIFEoper.INTforces.weights ;  % CECM weights

ndim = size(Xe,1) ; ngaus = length(weig) ; nnodeE = size(Xe,2)  ;
Ke = zeros(nnodeE*ndim,nnodeE*ndim) ;
nstrain = EIFEoper.MESH.nstrain ;
%dbstop('9')
for  g = 1:ngaus
    % Jacobian matrix at point g  (transformation of coordinates )
    % Matrix of derivatives for Gauss point "g" (polynomial shape functions)
    %  BeXi = dershapef(:,:,g) ;
    BeXi_transf = zeros(ndim,nnodeE) ;
    for idim = 1:ndim
        BeXi_transf(idim,:) =  EIFEoper.INTforces.Bgrad_transf{idim}(g,:) ;
    end
    % Jacobian Matrix
    Je = Xe*BeXi_transf' ;
    % JAcobian
    detJe = det(Je) ;
    % B-matrix for the parent element
    % --------------------------------
    BeXi =  EIFEoper.INTforces.Bgrad_cell{g} ;
    BeTILDE = inv(Je)'*BeXi ;
    Be_fine = QtransfB(BeTILDE,ndim) ;
    
    if nstrain == 4
        Be_fine = [Be_fine; zeros(1,size(Be_fine,2))] ;
    end
    Be = Be_fine*EIFEoper.INTforces.UdownsDEF_cell{g} ;    
    istrain = small2large(g,nstrain ) ;
    celas = EIFEoper.INTforces.MATPRO.celasglo(istrain,:) ;    
    Ke = Ke + weig(g)*detJe*(Be'*celas*Be) ;
end