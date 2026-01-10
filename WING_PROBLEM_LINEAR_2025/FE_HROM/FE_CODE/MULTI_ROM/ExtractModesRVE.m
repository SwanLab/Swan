function  DATAIN = ExtractModesRVE(DATA_TRAINING,DATAIN)
% Modal analysis domain-wise,
%dbstop('4')
if nargin ==0
    load('tmp.mat')
end
% if exist('ExtractDisplMatrix')~=2
%     addpath('FE_CODE') ;
% end
% if exist('BeamStiffnessMatrix')~=2
%     addpath('FE_CODE/BeamROM') ;
% end
if exist('SVDT')==0
    addpath('SVDlibrary')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT inputs
DATAIN= DefaultField(DATAIN,'NMODES_SHOW',[]) ;
DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;

DATAIN = DefaultField(DATAIN,'CUBATURE',[]) ;
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'CUBATURE',0) ; % Efficient integration scheme
DATAIN = DefaultField(DATAIN,'TOLERANCE_SVD_DISPLACEMENTS',[]) ;  % If not-empty, the number of modes are determined by this tolerance
%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%

% -------------------------------------
% CONSTRUCTING BASIS MATRICES
% ------------------------------------
DATAIN = DefaultField(DATAIN,'COMPUTE_MODES_AGAIN',1) ;
if DATAIN.COMPUTE_MODES_AGAIN == 1
    % ----------------------------------------------------------------------
    [BASES,DATA_REFMESH,~,~,stressDOM,MSG ]= BasisU_R_def_RVE(DATA_TRAINING,DATAIN) ;
    % -------------------------------------------------------------------------
    
    save(DATAIN.NAME_WS_MODES,'BASES','DATA_REFMESH','-append');
    
    DATAIN = DefaultField(DATAIN,'CUBATURE',[]); 
    DATAIN.CUBATURE =  DefaultField(DATAIN.CUBATURE,'INTERNAL_FORCES_WITHOUT_STRESS_MODES',0); 
    if DATAIN.CUBATURE.INTERNAL_FORCES_WITHOUT_STRESS_MODES == 1 % Internal force modes without stress modes
       DATAIN.NAMEWS_STRESS_SNAPSHOTS = [DATAIN.NAME_WS_MODES(1:end-4),'_stress.mat'] ; % New variable, 26-Jan-2020

        save(DATAIN.NAMEWS_STRESS_SNAPSHOTS,'stressDOM');
    else
        DATAIN.NAMEWS_STRESS_SNAPSHOTS = [] ; 
    end
    
else
    load(DATAIN.NAME_WS_MODES,'BASES','DATA_REFMESH') ;
    DATAIN = DefaultField(DATAIN,'NAMEWS_STRESS_SNAPSHOTS',[]);
    if ~isempty(DATAIN,'NAMEWS_STRESS_SNAPSHOTS')
        load(DATAIN.NAMEWS_STRESS_SNAPSHOTS,'stressDOM') ;
    end
end


if isempty(DATA_REFMESH.NODES_CORNERS)
    [DATAROM,BASES,MSG ]= RVEStiffnessMatrix(BASES,DATA_REFMESH,DATAIN,MSG) ;
else
    DATAROM = RVEStiffnessMatrix_PLATES(BASES,DATA_REFMESH,DATAIN) ;
end

save(DATAIN.NAME_WS_MODES,'DATAROM','-append');


%%% COMPUTING REDUCED (RVE STIFFNESS MATRIX) OF THE REFERENCE ELEMENT

DATAIN = DefaultField(DATAIN,'COMPUTE_STRESSES_AT_REDUCED_POINTS',1) ;

%if  DATAIN.COMPUTE_STRESSES_AT_REDUCED_POINTS == 1
%HROMVAR = ReducedSetIntegrationPoints(DATAIN,DATAROM,BASES,DATA_REFMESH) ;
[HROMVAR,MSG ]= ReducedSetIntegrationPoints(DATAIN,DATAROM,BASES,DATA_REFMESH,MSG) ;
% 
% % Reduced-order operators required for nonlinear analysis
% HROMVAR = ROMoperatNONL(DATAROM,BASES,DATA_REFMESH,DATAIN,HROMVAR) ; 

HROMVAR = ROMoperatNONL(DATAROM,BASES,DATA_REFMESH,DATAIN,HROMVAR) ; 


DATAROM.HROMVAR = HROMVAR ;
save(DATAIN.NAME_WS_MODES,'DATAROM','-append')
%end

% Nonlinear 

FILE_BATCH  = [cd,'/INFO.txt'] ;
fid =fopen(FILE_BATCH,'w');
for i = 1:length(MSG)
    fprintf(fid,[MSG{i},'\n']);
end
fod =fclose(fid);

open(FILE_BATCH);

