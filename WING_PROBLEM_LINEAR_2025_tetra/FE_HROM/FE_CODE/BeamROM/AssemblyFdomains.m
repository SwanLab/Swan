function [F,fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomains(DATAROM,MESH1D,DATAIN,FORCES_INPUT,ndim,DATA_REFMESH_glo)

if nargin == 0
    load('tmp.mat')
    %    DATAIN.INCLUDE_GRAVITY  = 1;
end



% T1/
DATAIN = DefaultField(DATAIN,'TRACTION_FORCES_DEFINED_BY_LINES',0);

VECTORIZED_CODE  =1;

if VECTORIZED_CODE == 0 | DATAIN.TRACTION_FORCES_DEFINED_BY_LINES == 0
    [F,fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= ...
        AssemblyFdomains_serial(DATAROM,MESH1D,DATAIN,FORCES_INPUT,ndim)  ;
else
    %  error('Something is wrong!')
    
    % We develop a new function for curved domains
    % --------------------------------------------
    DATA_REFMESH_glo{1} = DefaultField( DATA_REFMESH_glo{1},'RotationMatrixFace',[]) ;
    if isempty(DATA_REFMESH_glo{1}.RotationMatrixFace) || isempty(DATA_REFMESH_glo{1}.RotationMatrixFace{2})
        [F,fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= ...
            AssemblyFdomains_vector(DATAROM,MESH1D,DATAIN,FORCES_INPUT,ndim,DATA_REFMESH_glo)  ;
        
    else
        [F,fextBEAMr,rRB,fextDOMred,FORCES_2_PRINT]= ...
            AssemblyFdomains_vectorCURV(DATAROM,MESH1D,DATAIN,FORCES_INPUT,ndim,DATA_REFMESH_glo)  ;
    end
end
