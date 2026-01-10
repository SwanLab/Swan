function [F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    AssemblyFdomainsRVE_PLATE(DATAROM,MESH2D,DATAIN,FORCES_INPUT,ndim)

if nargin == 0
    load('tmp1.mat')
 %   DATAIN.INCLUDE_GRAVITY  = 1;
end



% T1/
DATAIN = DefaultField(DATAIN,'TRACTION_FORCES_DEFINED_BY_SURFACES',1);

VECTORIZED_CODE  =1;

if VECTORIZED_CODE == 0 | DATAIN.TRACTION_FORCES_DEFINED_BY_SURFACES == 0
    error('Option not implemented')
    %     [F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
    %     AssemblyFdomains_serial(DATAROM,MESH1D,DATAIN,FORCES_INPUT,ndim)  ;
else
    %  error('Something is wrong!')
    [F,fextRVEr,rRB,fextDOMred,FORCES_2_PRINT]= ...
        AssemblyFdomains_vectorRVE_plate(DATAROM,MESH2D,DATAIN,FORCES_INPUT,ndim)  ;
end


end

