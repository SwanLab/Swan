function [DISPLACEMENTS,REACTIONS,STRESSES] = ...
    DeterminationModes(DATAIN,dDOM,reactDOM,stressDOM,DATA_REFMESH,INDICES_PROJECTS)

if nargin == 0
    load('tmp4.mat')
end
DATAIN = DefaultField(DATAIN,'TypeMode','_BEAM') ; 

DATAIN = DefaultField(DATAIN,'nstepsGAUSS',[]) ;  
DATAIN = DefaultField(DATAIN,'nstepsNODES',[]) ;  

DATAIN.LEGEND_GRAPHS = ['Displacements',DATAIN.TypeMode] ;
nfigure.ERROR = 1;
nfigure.RIGHTSV = 1000;
LEGENDG =  ['Displacements',DATAIN.TypeMode] ;
COLOR = 'r' ;
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_DISPLACEMENTS(INDICES_PROJECTS) ;


DATAIN.PARTITION_COLUMNS = DATAIN.nstepsNODES ; 

DATAIN = DefaultField(DATAIN,'BeamProjectsLOC',[]) ; 
if ~isempty(DATAIN.BeamProjectsLOC)
    DATAIN.TOL_SVD_GLO = DATAIN.BeamProjectsLOC.TOL_SVD_GLO.DISP ; 
    ISBEAM = 1; 
else
    DATAIN.TOL_SVD_GLO = 0 ; % Relative, global error 
    ISBEAM = 0; 
end

[BasisUdef,SingVal_disp ]= ModesSVDplot(DATAIN,dDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;

DISPLACEMENTS.U = BasisUdef ;
DISPLACEMENTS.S = SingVal_disp ;


% REPEAT THE SAME PROCESS WITH REACTIONS
% -----------------------------------
if ISBEAM == 0
DATAIN.TOLERANCE_SVD_SLICE = cell(size(DATAIN.TOLERANCE_SVD_SLICE)) ; % We set to zero this tolerance, 
% because in principle we want all reaction and stress modes when dealing
% with nonbeam projects
end

DATAIN.LEGEND_GRAPHS = ['Reactions',DATAIN.TypeMode] ;
nfigure.RIGHTSV = 2000;
nfigure.ERROR = 1;
LEGENDG = ['Reactions',DATAIN.TypeMode] ;
COLOR = 'b' ;
DATAIN.PARTITION_COLUMNS = DATAIN.nstepsNODES ; 
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_REACTIONS(INDICES_PROJECTS) ;

[BasisRdef,SingVal_react ]= ModesSVDplot(DATAIN,reactDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;

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
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_STRESSES(INDICES_PROJECTS) ;
DATAIN.PLOT_MODES = 0 ;

[BasisS,SingVal_stress ]= ModesSVDplot(DATAIN,stressDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;


STRESSES.U = BasisS ;
STRESSES.S = SingVal_stress ;