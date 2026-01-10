function [DISPLACEMENTS,REACTIONS,STRESSES] = ...
    DeterminationModesRVE(DATAIN,dDOM,reactDOM,stressDOM,DATA_REFMESH,INDICES_PROJECTS)

if nargin == 0
    load('tmp.mat')
end
DATAIN = DefaultField(DATAIN,'TypeMode',' ESSENTIAL') ;
DATAIN.PLOT_RIGHT_SINGULAR_VECTORS = 0 ;
DATAIN.LEGEND_GRAPHS = ['Displacements',DATAIN.TypeMode] ;
nfigure.ERROR = 1;
nfigure.RIGHTSV = 1000;
LEGENDG =  ['Displacements',DATAIN.TypeMode] ;
COLOR = 'r' ;
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_DISPLACEMENTS(INDICES_PROJECTS) ;
DATAIN = DefaultField(DATAIN,'NMODES_trunc',[]) ; 
DATAIN.NMODES_trunc = DefaultField(DATAIN.NMODES_trunc,'DISPLACEMENTS',1e10) ; 
DATAIN.NMODES_TRUNCATE = DATAIN.NMODES_trunc.DISPLACEMENTS  ; 

% DATAIN.TOL_SVD_GLO.DISP = 1e-5 ; 
% DATAIN.TOL_SVD_GLO.REACTIONS = 1e-6 ; 
% DATAIN.TOL_SVD_GLO.STRESSES = 1e-6 ; 
DATAIN = DefaultField(DATAIN,'TOL_SVD_GLO',[]) ; 
DATAIN.TOL_SVD_GLO = DefaultField(DATAIN.TOL_SVD_GLO,'DISP',0) ; 
DATAIN.TOL_SVD_GLO = DefaultField(DATAIN.TOL_SVD_GLO,'REACTIONS',0) ; 
DATAIN.TOL_SVD_GLO = DefaultField(DATAIN.TOL_SVD_GLO,'STRESSES',0) ; 
TOL_SVD_GLO = DATAIN.TOL_SVD_GLO ; 

DATAIN.TOL_SVD_GLO = TOL_SVD_GLO.DISP ; 

[BasisUdef,SingVal_disp ]= ModesSVDplot(DATAIN,dDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;

PLOT_COVARIANCE_MATRIX = 0;
if PLOT_COVARIANCE_MATRIX == 1
    MATRIX = cell2mat(dDOM) ;
    COV =   PlotCovarianceMatrix(MATRIX) ;
end



DISPLACEMENTS.U = BasisUdef ;
DISPLACEMENTS.S = SingVal_disp ;


% REPEAT THE SAME PROCESS WITH REACTIONS
% ------------------

DATAIN.LEGEND_GRAPHS = ['Reactions',DATAIN.TypeMode] ;
nfigure.RIGHTSV = 2000;
nfigure.ERROR = 1;
LEGENDG = ['Reactions',DATAIN.TypeMode] ;
COLOR = 'b' ;
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_REACTIONS(INDICES_PROJECTS) ;
DATAIN.NMODES_trunc = DefaultField(DATAIN.NMODES_trunc,'REACTIONS',1e10) ; 
DATAIN.NMODES_TRUNCATE = DATAIN.NMODES_trunc.REACTIONS  ; 
DATAIN.TOL_SVD_GLO = TOL_SVD_GLO.REACTIONS ; 

[BasisRdef,SingVal_react ]= ModesSVDplot(DATAIN,reactDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;

REACTIONS.U = BasisRdef ;
REACTIONS.S = SingVal_react ;

if PLOT_COVARIANCE_MATRIX == 1
    MATRIX = cell2mat(reactDOM) ;
    COV =   PlotCovarianceMatrix(MATRIX) ;
end

% REPEAT THE SAME PROCESS WITH STRESSES
% ------------------

DATAIN.LEGEND_GRAPHS = ['Stresses',DATAIN.TypeMode] ;
nfigure.RIGHTSV = 3000;
nfigure.ERROR = 1;
LEGENDG = ['Stresses',DATAIN.TypeMode] ;;
COLOR = 'k' ;
DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_STRESSES(INDICES_PROJECTS) ;
DATAIN.PLOT_MODES = 0 ;
DATAIN.NMODES_trunc = DefaultField(DATAIN.NMODES_trunc,'STRESSES',1e10) ; 
DATAIN.NMODES_TRUNCATE = DATAIN.NMODES_trunc.STRESSES  ; 
DATAIN.TOL_SVD_GLO = TOL_SVD_GLO.STRESSES ; 

[BasisS,SingVal_stress ]= ModesSVDplot(DATAIN,stressDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS) ;


if PLOT_COVARIANCE_MATRIX == 1
    MATRIX = cell2mat(stressDOM) ;
    COV =   PlotCovarianceMatrix(MATRIX) ;
end

STRESSES.U = BasisS ;
STRESSES.S = SingVal_stress ;

end

function COV =  PlotCovarianceMatrix(SIGMA_SNAP,NFIG)
if nargin == 1
    NFIG = 28 ;
end

% Normalized ensemble
normSIGMA_SNAP = sqrt(sum((SIGMA_SNAP).*(SIGMA_SNAP)));
SIGMA_SNAPC = SIGMA_SNAP./repmat(normSIGMA_SNAP,size(SIGMA_SNAP,1),1);
marker_loc = '.';
labelMOD = 'NORM.';
labelTITLE = 'XSNAP = X_{\sigma}/norm(X_{\sigma})'  ;
COV = abs(SIGMA_SNAPC'*SIGMA_SNAPC) ;
% COV_COLOR = 256*COV;
% COV_COLOR_TRAC = (COV_TRAC-min_shear)*256/(1-min_shear);
figure(NFIG)

hold on
title('Colormap representation of the covariance matrix')
imagesc(COV)

end