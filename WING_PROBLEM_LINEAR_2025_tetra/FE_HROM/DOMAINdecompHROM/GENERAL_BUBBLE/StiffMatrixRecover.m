function Kstiff = StiffMatrixRecover(MATPRO,OPERFE) 



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
    celastST = ConvertBlockDiag(celastST) ; % Diagonal block matrix
    Kstiff = OPERFE.Bst'*(celastST*OPERFE.Bst);
else
    error('This option requires computing the stiffness matrix of the domain ')
end