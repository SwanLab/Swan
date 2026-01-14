function Kstiff = StiffMatrixRecoverLARGE_eifem(celastST,Bst,DATA,wSTs) 
 
if nargin == 0
    load('tmp3.mat')
end
 
    
    nF = size(celastST,2) ;
    for icomp = 1:nF
        icol = icomp:nF:size(celastST,1) ;
        celastST(icol,:) = bsxfun(@times,celastST(icol,:),wSTs) ;
    end
    
    
    ndim = DATA.MESH.ndim ; 
     FgradST = zeros(size(Bst,1),1) ;
 for idim = 1:ndim
     LOCROWS = idim:ndim^2:length(FgradST) ;
     FgradST(LOCROWS,:) = 1;
 end
    
    
    celastST = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;

    
    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
    Kstiff = Bst'*(celastST*Bst);
 