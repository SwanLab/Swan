function OFFLINE_STEPSbubble(TrainingTrajectoryINP,DATAoffline)
% Offline steps for constructing a hyperreduced-order  MULTISCALE MODEL
% Bubble modes option, see 
% 
% JAHO, 5-JUN-2023
%-------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
delete('OFFLINE.txt')
diary 'OFFLINE.txt'
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
% for controling SVD accuracy
DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
DATAoffline = DefaultField(DATAoffline,'errorSTRESS',0) ;
DATAoffline = DefaultField(DATAoffline,'errorFINT',1e-10) ;
DATAoffline = DefaultField(DATAoffline,'errorECM',0) ;
DATAoffline = DefaultField(DATAoffline,'errorPK2stress_basis',DATAoffline.errorSTRESS) ;
DATAoffline = DefaultField(DATAoffline,'errorREACTIVE_FORCES',0) ;

%
% READING INPUT DATA
DATAcommon = feval(TrainingTrajectoryINP) ; 
% Offline steps for constructing a hyperreduced-order  MULTISCALE MODEL
% JAHO, 30-JAN-2023
%-------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
delete('OFFLINE.txt')
diary 'OFFLINE.txt'
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ;
end
% for controling SVD accuracy
DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
DATAoffline = DefaultField(DATAoffline,'errorSTRESS',0) ;
DATAoffline = DefaultField(DATAoffline,'errorFINT',1e-10) ;
DATAoffline = DefaultField(DATAoffline,'errorECM',0) ;
DATAoffline = DefaultField(DATAoffline,'errorPK2stress_basis',DATAoffline.errorSTRESS) ;
DATAoffline = DefaultField(DATAoffline,'errorREACTIVE_FORCES',0) ;

%
% READING INPUT DATA
DATAcommon = feval(TrainingTrajectoryINP) ;
NAME_BASE =[DATAcommon.NameParamStudyLOC,'_param_'];
NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];

% -------------------------------------------
% 1) RETRIEVING DISPLACEMENT MODES AND MESH INFORMATION
% ------------------------------------------
[INFO_RVE,MESHdom,SNAPdisp,BasisU,DATA,OTHER_output,OPERFEaux,MATPRO,Fbody,Ftrac] ...
    = GetDisplacAndMesh_1domBUB(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline) ;

% -------------------------------------------
% 2) RETRIEVING STRESSES AND REACTIVE FORCES
% ------------------------------------------
[SNAPreact,BasisStwo,BasisPone,DATA.NAME_BASE]  = GetStressesAndReactForces_1domBUB(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
    INFO_RVE,MESHdom,BasisU,DATA,OPERFEaux,MATPRO,Fbody,Ftrac,OTHER_output) ;

OPERFEaux = [] ;


% -------------------------------------------------------------------------------------------------------------------------------------------------
% 3) FE OPERATORS AND OTHER VARIABLES RELATED WITH THE FE DISCRETIZATION OF THE DOMAIN (SUCH AS GLOBAL SHAPE FUNCTIONS, Global B-matrix operator)
% --------------------------------------------------------------------------------------------------------------------------------------------------
[OPERFE,MESHdom,DATA] = FEoperators_1dom(DATA,MESHdom,MATPRO)   ;

% 4) MODES

DATAoffline= DefaultField(DATAoffline,'ReactiveModesFromInterfaceDisplacements',0) ;


if DATAoffline.ReactiveModesFromInterfaceDisplacements == 0
    % VERSION BEFORE 1-APRIL-2023
    % IN THIS VERSION, REACTIVE MODES ARE CALCULATED FIRST (TRUNCATED SVD)
    [PhiDEF,PsiDEFf,Vall,Mintf,Mdom,PhiRB,PsiRBf,MESHdom,DATA]  =...
        MODESandROMoper_1dom(OPERFE,MESHdom,DATA,SNAPreact,BasisStwo,...
        SNAPdisp,DATAcommon,DATAoffline) ;
elseif DATAoffline.ReactiveModesFromInterfaceDisplacements == 1
    % TENTATIVE VERSION 1-APRIL-2023
    % WE FOUND IN
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx
    % that the version used in giuliudory2023multiscale.pdf is prone to
    % instabilities. The version implemented here attempts to eliminate
    % such instabilities
    [PhiDEF,PsiDEFf,Vall,Mintf,Mintf_INV,Mdom,PhiRB,PsiRBf,MESHdom,DATA]  =...
        MODESandROMoper_EIFE(OPERFE,MESHdom,DATA,SNAPreact,BasisStwo,...
        SNAPdisp,DATAcommon,DATAoffline) ;
    
else
    error('Option not implemented')
end




% 5) HYPERREDUCTION (SELECTION POINTS AND WEIGHTS)
%-------------------

%----------------------------------
DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',[]) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'INTERNAL_FORCES',0) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'MASS_MATRIX',0) ;

CECM= ContinuousECM_multi(OPERFE,PhiDEF,BasisPone,DATA,DATAcommon,MESHdom,DATAoffline,MATPRO,PhiRB) ;

MODES.PhiDEF = PhiDEF;
MODES.PsiDEFf = PsiDEFf;
MODES.Vall = Vall;
MODES.PhiRB = PhiRB;
MODES.PsiRBf = PsiRBf;
MODES.BasisPone = BasisPone;
MODES.BasisStwo = BasisStwo;


GEOMATRICES.Mintf = Mintf;
%GEOMATRICES.MdomFF = MdomFF;
GEOMATRICES.Mdom = Mdom;

DATA_misc = DATA;



save(NAMEOFFLINEstore,'MODES','GEOMATRICES','MESHdom','DATA_misc','DATAoffline','MATPRO','OPERFE','CECM','DATAcommon')



diary off

% 
% 
% NAME_BASE =[DATAcommon.NameParamStudyLOC,'_param_'];
% NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
% NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];
% 
% % -------------------------------------------
% % 1) RETRIEVING DISPLACEMENT MODES AND MESH INFORMATION
% % ------------------------------------------
% [INFO_RVE,MESHdom,SNAPdisp,BasisU,DATA,OTHER_output,OPERFEaux,MATPRO,Fbody,Ftrac] ...
%     = GetDisplacAndMesh_1dom(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline) ;
% 
% % -------------------------------------------
% % 2) RETRIEVING STRESSES AND REACTIVE FORCES
% % ------------------------------------------
% [SNAPreact,BasisStwo,BasisPone,DATA.NAME_BASE]  = GetStressesAndReactForces_1dom(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
%     INFO_RVE,MESHdom,BasisU,DATA,OPERFEaux,MATPRO,Fbody,Ftrac,OTHER_output) ;
% 
% OPERFEaux = [] ;
% 
% 
% % -------------------------------------------------------------------------------------------------------------------------------------------------
% % 3) FE OPERATORS AND OTHER VARIABLES RELATED WITH THE FE DISCRETIZATION OF THE DOMAIN (SUCH AS GLOBAL SHAPE FUNCTIONS, Global B-matrix operator)
% % --------------------------------------------------------------------------------------------------------------------------------------------------
% [OPERFE,MESHdom,DATA] = FEoperators_1dom(DATA,MESHdom,MATPRO)   ;
% 
% % 4) MODES
% 
% DATAoffline= DefaultField(DATAoffline,'ReactiveModesFromInterfaceDisplacements',0) ;
% 
% 
% if DATAoffline.ReactiveModesFromInterfaceDisplacements == 0
%     % VERSION BEFORE 1-APRIL-2023
%     % IN THIS VERSION, REACTIVE MODES ARE CALCULATED FIRST (TRUNCATED SVD)
%     [PhiDEF,PsiDEFf,Vall,Mintf,Mdom,PhiRB,PsiRBf,MESHdom,DATA]  =...
%         MODESandROMoper_1dom(OPERFE,MESHdom,DATA,SNAPreact,BasisStwo,...
%         SNAPdisp,DATAcommon,DATAoffline) ;
% elseif DATAoffline.ReactiveModesFromInterfaceDisplacements == 1
%     % TENTATIVE VERSION 1-APRIL-2023
%     % WE FOUND IN
%     % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx
%     % that the version used in giuliudory2023multiscale.pdf is prone to
%     % instabilities. The version implemented here attempts to eliminate
%     % such instabilities
%     [PhiDEF,PsiDEFf,Vall,Mintf,Mintf_INV,Mdom,PhiRB,PsiRBf,MESHdom,DATA]  =...
%         MODESandROMoper_EIFE(OPERFE,MESHdom,DATA,SNAPreact,BasisStwo,...
%         SNAPdisp,DATAcommon,DATAoffline) ;
%     
% else
%     error('Option not implemented')
% end
% 
% 
% 
% 
% % 5) HYPERREDUCTION (SELECTION POINTS AND WEIGHTS)
% %-------------------
% 
% %----------------------------------
% DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',[]) ;
% DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'INTERNAL_FORCES',0) ;
% DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'MASS_MATRIX',0) ;
% 
% CECM= ContinuousECM_multi(OPERFE,PhiDEF,BasisPone,DATA,DATAcommon,MESHdom,DATAoffline,MATPRO,PhiRB) ;
% 
% MODES.PhiDEF = PhiDEF;
% MODES.PsiDEFf = PsiDEFf;
% MODES.Vall = Vall;
% MODES.PhiRB = PhiRB;
% MODES.PsiRBf = PsiRBf;
% MODES.BasisPone = BasisPone;
% MODES.BasisStwo = BasisStwo;
% 
% 
% GEOMATRICES.Mintf = Mintf;
% %GEOMATRICES.MdomFF = MdomFF;
% GEOMATRICES.Mdom = Mdom;
% 
% DATA_misc = DATA;
% 
% 
% 
% save(NAMEOFFLINEstore,'MODES','GEOMATRICES','MESHdom','DATA_misc','DATAoffline','MATPRO','OPERFE','CECM','DATAcommon')
% 
% 
% 
% diary off
% 
% 
