function [DISPLACEMENTS,REACTIONS,STRESSES,MSG] = ...
    DeterminationModesUnified(DATAIN,dDOM,reactDOM,stressDOM,DATA_REFMESH,MSG)
% Copy of   DeterminationModes. New approach (10-Dec-2018)
% 

if nargin == 0
    load('tmp2.mat')
end
 
DATAIN = DefaultField(DATAIN,'TypeMode','') ; 

DATAIN = DefaultField(DATAIN,'nstepsGAUSS',[]) ;  
DATAIN = DefaultField(DATAIN,'nstepsNODES',[]) ;  

DATAIN.LEGEND_GRAPHS = ['Displacements',DATAIN.TypeMode] ;
nfigure.ERROR = 1;
nfigure.RIGHTSV = 1000;
LEGENDG =  ['Displacements',DATAIN.TypeMode] ;
COLOR = 'r' ;
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_DISPLACEMENTS ;


DATAIN.PARTITION_COLUMNS = DATAIN.nstepsNODES ; 

DATAIN = DefaultField(DATAIN,'BeamProjectsLOC',[]) ; 

DATAIN.TOL_SVD_GLO = 0 ; % Relative, global error 
%ISBEAM = 1; 


%%% DISPLACEMENTS
% ---------------------
if isfield(DATAIN,'TOLERANCE_SVD_SLICE')
    % Beam problems
    RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_SLICE ;
elseif  isfield(DATAIN,'TOLERANCE_SVD_RVES_DISP')
    % RVEs problems 
     RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_RVES_DISP ;
else
    RELATIVE_TOLERANCE = cell(size(dDOM)) ; 
    RELATIVE_TOLERANCE(:) = {0} ; 
end
 
 MSG{end+1} = '------------------------------------------'  ; 
 MSG{end+1} = 'SVD modes displacements (before alignment)'  ; 
  MSG{end+1} = '------------------------------------------'  ; 

 
[BasisUdef,SingVal_disp,MSG ]= ModesSVDplotUnified(DATAIN,dDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,RELATIVE_TOLERANCE,MSG) ;

DISPLACEMENTS.U = BasisUdef ;
DISPLACEMENTS.S = SingVal_disp ;


% REPEAT THE SAME PROCESS WITH REACTIONS
% -----------------------------------
 
DATAIN.LEGEND_GRAPHS = ['Reactions',DATAIN.TypeMode] ;
nfigure.RIGHTSV = 2000;
nfigure.ERROR = 1;
LEGENDG = ['Reactions',DATAIN.TypeMode] ;
COLOR = 'b' ;
DATAIN.PARTITION_COLUMNS = DATAIN.nstepsNODES ; 
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_REACTIONS ; %(INDICES_PROJECTS) ;

%[BasisRdef,SingVal_react ]= ModesSVDplot(DATAIN,reactDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;

%RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_SLICE_REACT ; 
%%% REACTIONS
% ---------------------
if isfield(DATAIN,'TOLERANCE_SVD_SLICE_REACT')
    % Beam problems
    RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_SLICE_REACT ;
elseif  isfield(DATAIN,'TOLERANCE_SVD_RVES_REACT')
    % RVEs problems 
     RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_RVES_REACT ;
else
    RELATIVE_TOLERANCE = cell(size(dDOM)) ; 
    RELATIVE_TOLERANCE(:) = {0} ; 
 end
 MSG{end+1} = '------------------------------------------'  ; 
 MSG{end+1} = 'SVD modes reactions  '  ; 
  MSG{end+1} = '------------------------------------------'  ; 
[BasisRdef,SingVal_react,MSG ]= ModesSVDplotUnified(DATAIN,reactDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,...
    RELATIVE_TOLERANCE,MSG) ;



REACTIONS.U = BasisRdef ;
REACTIONS.S = SingVal_react ;

% REPEAT THE SAME PROCESS WITH STRESSES
% ------------------

DATAIN.LEGEND_GRAPHS = ['Stresses',DATAIN.TypeMode] ;
nfigure.RIGHTSV = 3000;
nfigure.ERROR = 1;
LEGENDG = ['Stresses',DATAIN.TypeMode] ;;
COLOR = 'k' ;
DATAIN.PARTITION_COLUMNS = DATAIN.nstepsGAUSS ; 
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_STRESSES ; %(INDICES_PROJECTS) ;
DATAIN.PLOT_MODES = 1 ;
DATAIN.ISGAUSS = 1 ; 


%[BasisS,SingVal_stress ]= ModesSVDplot(DATAIN,stressDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;

%RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_SLICE_STRESS ;

%%% STRESSES
% ---------------------
if isfield(DATAIN,'TOLERANCE_SVD_SLICE_STRESS')
    % Beam problems
    RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_SLICE_STRESS ;
elseif  isfield(DATAIN,'TOLERANCE_SVD_RVES_STRESS')
    % RVEs problems 
     RELATIVE_TOLERANCE = DATAIN.TOLERANCE_SVD_RVES_STRESS ;
else
    RELATIVE_TOLERANCE = cell(size(dDOM)) ; 
    RELATIVE_TOLERANCE(:) = {0} ; 
 end

 MSG{end+1} = '------------------------------------------'  ; 
 MSG{end+1} = 'SVD modes stresses '  ; 
  MSG{end+1} = '------------------------------------------'  ; 
[BasisS,SingVal_stress,MSG ]= ModesSVDplotUnified(DATAIN,stressDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,...
    RELATIVE_TOLERANCE,MSG) ;



STRESSES.U = BasisS ;
STRESSES.S = SingVal_stress ;

DATAIN = DefaultField(DATAIN,'UseSnapshotStressesInComputingBasisF',0) ; %  % Variable created 17-Apr-2020. If = 1, then BasisF is computed with 
% the original stress snapshots.  Therefore, stressDOM is to be stored in
% memory without any alteration whatsoever 
if DATAIN.UseSnapshotStressesInComputingBasisF == 1
    STRESSES.SNAPSHOTS = cell2mat(stressDOM) ; 
else
    STRESSES.SNAPSHOTS  = [] , 
end
