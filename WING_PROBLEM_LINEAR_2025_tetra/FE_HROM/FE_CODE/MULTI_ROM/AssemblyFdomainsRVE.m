function [F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomainsRVE(DATAROM,MESH2D,DATAIN,FORCES_INPUT,ndim,DOFsKEEP,DATA_REFMESH_glo)

if nargin == 0
    load('tmp2.mat')
    %    DATAIN.INCLUDE_GRAVITY  = 1;
end



% T1/
DATAIN = DefaultField(DATAIN,'TRACTION_FORCES_DEFINED_BY_SURFACES',1);
DATAIN.nRB = size(DATAROM{1}.BasisIntRB{1},2) ;

VECTORIZED_CODE  =1;

if VECTORIZED_CODE == 0 | DATAIN.TRACTION_FORCES_DEFINED_BY_SURFACES == 0
    error('Option not implemented')
    %     [F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    %     AssemblyFdomains_serial(DATAROM,MESH1D,DATAIN,FORCES_INPUT,ndim)  ;
else
    %  error('Something is wrong!')
    % FORCES_INPUT = DefaultField(FORCES_INPUT,'LOAD',[]) ;
    % if ~isempty(FORCES_INPUT.LOAD) || DATAIN.INCLUDE_GRAVITY == 1
    
    DATA_REFMESH_glo{1} = DefaultField( DATA_REFMESH_glo{1},'RotationMatrixFace',[]) ;
    
    ROTATIONS = DATA_REFMESH_glo{1}.RotationMatrixFace ;
    ISROTATED = 0 ;
    if ~isempty(ROTATIONS)
        [AAA, BBB] = cellfun(@size,ROTATIONS) ;
        if any(AAA)
            ISROTATED = 1;
        end
    end
    
    if ISROTATED == 0
        
        [F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
            AssemblyFdomains_vectorRVE(DATAROM,MESH2D,DATAIN,FORCES_INPUT,ndim,DOFsKEEP)  ;
    else
        [F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
            AssemblyFdomains_vectorRVE_curved(DATAROM,MESH2D,DATAIN,FORCES_INPUT,ndim,DOFsKEEP,ROTATIONS)  ;
    end
    
end


end


