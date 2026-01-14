function [MATPRO] = SmallStrainElasticityPROP(MESH,typePROBLEM,PROPMAT)

ndim = size(MESH.COOR,2)  ;

MESH = DefaultField(MESH,'nstrain',[]) ;
if isempty(MESH.nstrain)
    if ndim==2
        nstrain = 3;
    else
        nstrain = 6 ;
        typePROBLEM ='3D' ;
    end
else
    nstrain = MESH.nstrain ;
end

nelem = size(MESH.MaterialType,1) ;
MATPRO.celasglo = zeros(nstrain,nstrain,nelem) ;  % Global array of elasticity matrices
MATPRO.dens = zeros(nelem,1) ;
%celasgloINV = zeros(6,6,nelem) ;
for imat = 1:length(PROPMAT)
    celas3D =PROPMAT(imat).ElasticityMatrix ; %
    INVcelas3D = inv(celas3D) ;
    ELEMS = find(MESH.MaterialType == imat) ;
    switch typePROBLEM
        case 'pstrain'
            if nstrain == 3
                rowcol = [1 2 6] ;
            elseif nstrain == 4
                rowcol = [1 2 6 3] ;
            else
                error('Option not implemented')
            end
            celas = celas3D(rowcol,rowcol) ;
            
        case 'pstress'
            if nstrain == 3
                rowcol = [1 2 6] ;
                celasINV3D = inv(celas3D) ;
                celasINV = celasINV3D(rowcol,rowcol) ;
                celas = inv(celasINV) ;
            else
                error('Option not implemented')
            end
        case '3D'
            celas = celas3D ;
    end
    for eLOC=1:length(ELEMS)
        e = ELEMS(eLOC) ;
        MATPRO.celasglo(:,:,e) = celas ;
        MATPRO.dens(e) = PROPMAT(imat).Density ;
        
    end
end
