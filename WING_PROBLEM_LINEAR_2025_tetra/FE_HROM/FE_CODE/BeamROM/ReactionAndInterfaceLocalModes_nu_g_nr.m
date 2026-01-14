function [BasisRdef,BasisINT,MSG,BasisUdef,SingVal_Udef] = ...
    ReactionAndInterfaceLocalModes_nu_g_nr(BasisUdef,BasisRdef,f1,f2,TOL_SINGULAR_VALUES_Hqr,...
    nBASES_BEAM,DATA_REFMESH,Vrb,M,DATAIN,SingVal_Udef,SingVal_Rdef,BasisUrb,Mdom,MSG)
if nargin == 0
    load('tmp.mat')
end

% copy of ReactionAndInterfaceLocalModes_UNIFIED.m  (Modification 13-May-2020 to allow the possibility of
% having more displacement modes than reaction modes)

f = [f1;f2] ;

% Candidates to be interface displacement modes
% ----------------------------------------------
% Rotation
% ---------
if ~isempty(DATA_REFMESH.RotationMatrixFace)
    iface = 2;
    DATAIN.ROTATION_LOC = DATA_REFMESH.RotationMatrixFace{iface} ;
    R = DATAIN.ROTATION_LOC ;
else
    R = [] ;
end


% Selection of reaction modes
% ---------------------------
nrigidB = size(Vrb,2) ;

% Determination of aligned displacement method
% Modified  on May-13th-2020, to filter out modes that  do not
% contribute to the work done by the  unit cell
DATAIN.FILTER_ORTHOGONAL_REACTION_DISPLACEMENT_MODES = 1 ;
DATAIN = DefaultField(DATAIN,'ANGLE_FILTER_ORTHOGONAL_REACTION_DISPLACEMENT_MODES',89) ; % Degrees
TOL_SVD = 0 ;
ANGLE = DATAIN.ANGLE_FILTER_ORTHOGONAL_REACTION_DISPLACEMENT_MODES ;
[RA,RB,ANGLES ]= AlignedSubspaces(BasisUdef,BasisRdef,TOL_SVD,ANGLE) ;
BasisUdef_eff = RA;
BasisRdef = RB ;
SingVal_Udef = ones(size(RA,2),1) ;
SingVal_Rdef = ones(size(RB,2),1) ;

 % NEW CONCEPT --> Effective straining modes --> BasisUdef_eff
 
MSG{end+1}  ='--------------------------------------' ; 
MSG{end+1}  = ['NUMBER OF DISPLACEMENT (strain) MODES = ',num2str(size(BasisUdef,2))];
MSG{end+1}  = ['NUMBER OF reactions = effective (strain) MODES = ',num2str(size(BasisUdef_eff,2))];
MSG{end+1}  = ['ALIGNMENT ANGLES = ',num2str(ANGLES(1:size(RB,2))')];



% PLOT AGAING DISPL. MODES

BasisUdefPlot = BasisUdef_eff ;

% end

CNref = DATA_REFMESH.CN  ;
COOR = DATA_REFMESH.COOR  ;
NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),'DISP_ALIGNED' ];
DATA= [] ;
LEGENDG = 'DISP.' ;

MSG{end+1} = '-----------------------------------';
MSG{end+1} = 'DISPLACEMENT MODES AFTER ALIGNMENT';
MSG{end+1} = '-----------------------------------';
MSG =   GidPostProcessModes_dom(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdefPlot,DATA_REFMESH.posgp,...
    NAME_MODES,DATA,LEGENDG,MSG);
MSG{end+1} = '-----------------------------------';


[Vdef_Morth,Vrb_Morth] = CandidatesInterfaceModes_UNIFIED(BasisUdef_eff,f1,f2,SingVal_Udef,...
    DATAIN,M,Vrb,BasisRdef,BasisUrb) ;


DATAIN = DefaultField(DATAIN,'TOL_DeformationalInterfaceModes_AlignmentMethod',1e-3) ;
DATAIN = DefaultField(DATAIN,'KINEMATIC_CONSTRAINTS_MODES',[]) ;
DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'DEIM_BASED_METHOD',1) ;

% warning('This implementation is not reliable')
DATAIN.RotationMatrixLocal{1} = [] ;
DATAIN.RotationMatrixLocal{2} = R' ;
iface = 1;
[BasisINT,~,MSG ]=  DeformModesInterface_AlignmentMethod(BasisRdef,f1,f2,...
    Vrb,M,DATAIN,BasisUdef_eff,SingVal_Udef,SingVal_Rdef,iface,MSG,Vdef_Morth) ;
% TEXTB ={} ;



