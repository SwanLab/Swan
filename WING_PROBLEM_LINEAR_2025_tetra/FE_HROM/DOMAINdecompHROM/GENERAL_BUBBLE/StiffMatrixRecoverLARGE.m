function Kstiff = StiffMatrixRecoverLARGE(MATPRO,OPERFE,MESH,DATA) 
 
if nargin == 0
    load('tmp3.mat')
end


disp(['Recovering stiffness matrix...'])

if isfield(MATPRO,'Celas')
    MATPRO.celasglo = MATPRO.Celas ;
end


if ~isempty(MATPRO.celasglo)
    %     disp('Checking accuracy CECM points ')
    %     disp('-------------------------------------------')
    %   disp('Coarse-scale stiffness matrix ')
    celastST = MATPRO.celasglo ;
    
     
    
    nF = size(MATPRO.celasglo,2) ;
    for icomp = 1:nF
        icol = icomp:nF:size(celastST,1) ;
        celastST(icol,:) = bsxfun(@times,celastST(icol,:),OPERFE.wSTs) ;
    end
    
    
    ndim = DATA.MESH.ndim ; 
     FgradST = zeros(size(OPERFE.Bst,1),1) ;
 for idim = 1:ndim
     LOCROWS = idim:ndim^2:length(FgradST) ;
     FgradST(LOCROWS,:) = 1;
 end
    
    
    celastST = CelasLARGEmat_allgauss(celastST,FgradST,ndim) ;

    
    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
    Kstiff = OPERFE.Bst'*(celastST*OPERFE.Bst);
else
    error('This option requires computing the stiffness matrix of the domain ')
end