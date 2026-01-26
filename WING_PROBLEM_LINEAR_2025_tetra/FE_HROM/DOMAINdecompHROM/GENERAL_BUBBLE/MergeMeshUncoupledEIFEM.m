function MESH = MergeMeshUncoupledEIFEM(MESH,DATAcommon)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/06_UNCOUPLED.mlx
% JAHO, 23-Apr-2024, Balmes 185, Barcelona
% ------------------------------------------------------------
if nargin == 0
    load('tmp2.mat')
end

% INPUT DATA
% MESH.CN
% MESH.COOR
DATAcommon = DefaultField(DATAcommon,'UNCOUPLED_EIFEM',[]) ;
DATAcommon.UNCOUPLED_EIFEM = DefaultField(DATAcommon.UNCOUPLED_EIFEM,'INDICES_TYPEstructure',{[1],[2]}) ;

INDT = DATAcommon.UNCOUPLED_EIFEM.INDICES_TYPEstructure ;
if length(DATAcommon.UNCOUPLED_EIFEM.INDICES_TYPEstructure) == 2
    % 2D PROBLEMS
    % INITIAL VERSION
    if length(INDT{1}) >1
        error('option not implemented yet: only one type of structure per direction')
        
    else
        MESH = MergeMeshUncoupledEIFEM_2D_1mat(MESH,INDT,DATAcommon);
    end
    
    
    
    
elseif length(DATAcommon.UNCOUPLED_EIFEM.INDICES_TYPEstructure) == 2
    % 3D PROBLEMS
    error('Option not implemented yet')
end


