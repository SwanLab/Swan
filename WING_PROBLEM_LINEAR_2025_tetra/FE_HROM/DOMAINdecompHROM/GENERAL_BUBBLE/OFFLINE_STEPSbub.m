function OFFLINE_STEPSbub(TrainingTrajectoryINP,DATAoffline)
% Offline steps for constructing a hyperreduced-order  MULTISCALE MODEL
% Bubble modes
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/05_Assessment.mlx
% JAHO, 30-JAN-2023/23-Sept-2023/17-Oct-2023
%-------------------------------------------------------------
if nargin == 0
    load('tmp.mat')
end
delete('OFFLINE.txt')
diary 'OFFLINE.txt'
if ~exist('Quadratic3N')
    addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ; %
end
% for controling SVD accuracy
% DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
% DATAoffline = DefaultField(DATAoffline,'errorSTRESS',0) ;
DATAoffline = DefaultField(DATAoffline,'errorFINT',1e-10) ;
% DATAoffline = DefaultField(DATAoffline,'errorECM',0) ;
% DATAoffline = DefaultField(DATAoffline,'errorPK2stress_basis',DATAoffline.errorSTRESS) ;
% DATAoffline = DefaultField(DATAoffline,'errorREACTIVE_FORCES',0) ;
% DATAoffline = DefaultField(DATAoffline,'errorSVD_SELF_EQUILIBRATED',1e-6) ;

DATAoffline = DefaultField(DATAoffline,'TOLSVD_complementary_modes',1e-3) ;



DATAoffline = DefaultField(DATAoffline,'TOLSVD_complementary_modes_stresses',1e-10) ;
DATAoffline = DefaultField(DATAoffline,'BASIS_STRESSES_FROM_PROJECTED_DISPLACEMENTS',1) ;


% DATAoffline.TOLSVD_complementary_modes = 1e-3 ; % displacements
% DATAoffline.TOLSVD_complementary_modes_stresses = 1e-3 ;
% DATAoffline.BASIS_STRESSES_FROM_PROJECTED_DISPLACEMENTS = 0 ;



%
% READING INPUT DATA
DATAcommon = feval(TrainingTrajectoryINP) ;
DATAcommon = DefaultField(DATAcommon,'NameFileMeshDATA',[] ) ;
DATAoffline = DefaultField(DATAoffline,'NAMEMESH', DATAcommon.NameFileMeshDATA);
DATAcommon.NameFileMeshDATA = DATAoffline.NAMEMESH ;

NAME_BASE =[DATAcommon.NameParamStudyLOC,'_param_'];
NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];

% -------------------------------------------
% 1) RETRIEVING DISPLACEMENT MODES AND MESH INFORMATION
% ------------------------------------------
DATAoffline.ReturnFullDisplacementMatrix = 1;
DATAoffline.USE_BASIS_MATRIX_DISPLACEMENTS_AS_SNAPSHOTS = 2;
[INFO_RVE,MESHdom,SNAPdisp,DATA,OTHER_output,OPERFEaux,MATPRO,Fbody,Ftrac,SNAPdisp_AUX,DATAoffline] ...
    = GetDisplacAndMesh_sevDOM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline) ;


% -------------------------------------------------------------------------------------------------------------------------------------------------
% 2) FE OPERATORS AND OTHER VARIABLES RELATED WITH THE FE DISCRETIZATION OF THE DOMAIN (SUCH AS GLOBAL SHAPE FUNCTIONS, Global B-matrix operator)
% --------------------------------------------------------------------------------------------------------------------------------------------------
[OPERFE,MESHdom,DATA] = FEoperators_1dom(DATA,MESHdom,MATPRO)   ;


% DEFORMATIONAL MODES (PURGE RIGID-BODY COMPONENT).
% -------------------

DATAoffline = DefaultField(DATAoffline,'METHOD_INVARIANT_DEF_MODES',[]);

%DATAoffline.METHOD_INVARIANT_DEF_MODES = DefaultField(DATAoffline.METHOD_INVARIANT_DEF_MODES,'STAGE_TRAINING',1);

 [BasisUdeform,Mdom,PhiRB,MESHdom,DATA,PsiRBf,Mintf,DATAoffline]  =...
        DeformModesNONL_general(OPERFE,MESHdom,DATA,...
        SNAPdisp,DATAcommon,DATAoffline,SNAPdisp_AUX) ;
    

% -------------------------------------------
% 2) RETRIEVING STRESSES AND REACTIVE FORCES  (project-wise)
% ------------------------------------------
DATAoffline = DefaultField(DATAoffline,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; % 6-Jan-2023,  see
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx

if DATAoffline.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
    AsnapSTRESSinel = [] ;
    [AsnapREAC ,AsnapSTRESS ,~,NAME_BASE]  = GetStressesAndReactForces_bub(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
        INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
else
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
    % Separate treatment nonlinear stresses
    [AsnapREAC ,AsnapSTRESS ,~,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel]  = GetStressesAndReactForces_bubNECM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
        INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
    
end




% 3) Self-equilibrated modes and Bubble Modes
[PsiSEf,PhiDEF,GammaBUB] =     BubbleModes_BasicSE(PsiRBf,Mintf,AsnapREAC,DATA,DATAoffline,MESHdom,BasisUdeform,Mdom)  ;

DATAoffline =DefaultField(DATAoffline,'CHECK_ACCURACY_PROJECTED_STRESSES_AFTER_BUBBLE_MODES',0) ;

if DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES_AFTER_BUBBLE_MODES == 1
    
    BasisUdeform_bubble.PhiDEFbs = PhiDEF ;
    BasisUdeform_bubble.PhiDEFcomp = GammaBUB;
    
    
    DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES = 1 ;
    [~ ,~ ,~,~]  = GetStressesAndReactForces_bub(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
        INFO_RVE,BasisUdeform_bubble,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
    
end



% 4) Stress modes
% ---------------------------------------------------------------------
DATAoffline = DefaultField(DATAoffline,'StressModesComputation_viaCOEFFICIENTS_SVD_DEFORMATION',0) ;
if DATAoffline.StressModesComputation_viaCOEFFICIENTS_SVD_DEFORMATION == 1
    % Seldom used. I dare to say it may be deleted (15-03-2024)
    [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_bub(AsnapSTRESS,DATA,DATAoffline,MESHdom,OPERFE,BasisUdeform)  ;
else
    % THIS IS THE STANDARD OPTION (7-JAN-2024), FORGET ABOUT THE ABOVE ONE
    if isempty(AsnapSTRESSinel)
        [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_viaSVD(AsnapSTRESS,DATA,DATAoffline,MESHdom,OPERFE,BasisUdeform)  ;
    else
        % 7-Jan-2024, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        % CECM via inelastic stresses
        [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_viaSVDinel(AsnapSTRESSinel,DATA,DATAoffline,MESHdom,OPERFE,BasisUdeform,AsnapSTRESSel)  ;
    end
end



% -------------------------------
% 5) Interface displacements
% --------------------------------
%GammaBUB = BasisUdeform.BUBBLE ;
% = BasisUdeform.BASIC;
% PsiDEFf = PsiDEFf.BASIC;

[Vall,MESHdom,DATA] = FictInterfaceDISP(PhiDEF,PsiSEf,MESHdom,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline);


%
% DATAoffline= DefaultField(DATAoffline,'ReactiveModesFromInterfaceDisplacements',0) ;
% DATAoffline= DefaultField(DATAoffline,'MethodDeformationalModesDirectly_OnlyElasticRange',0) ;
%
% if  DATAoffline.MethodDeformationalModesDirectly_OnlyElasticRange ==1
%     % 23-Aug-2023, see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/103_EIFEM_PLATES/02_CONV_SIZE_TRAC.mlx
%      [BasisUdeform,PsiDEFf,Vall,Mintf,Mdom,PhiRB,PsiRBf,MESHdom,DATA]  =...
%             MODESandROMoper_elastic(OPERFE,MESHdom,DATA,SNAPreact,BasisStwo,...
%             SNAPdisp,DATAcommon,DATAoffline,MATPRO) ;
%
% else
%
%     if DATAoffline.ReactiveModesFromInterfaceDisplacements == 0
%         % VERSION BEFORE 1-APRIL-2023
%         % IN THIS VERSION, REACTIVE MODES ARE CALCULATED FIRST (TRUNCATED SVD)
%         [BasisUdeform,PsiDEFf,Vall,Mintf,Mdom,PhiRB,PsiRBf,MESHdom,DATA]  =...
%             MODESandROMoper_1dom(OPERFE,MESHdom,DATA,SNAPreact,BasisStwo,...
%             SNAPdisp,DATAcommon,DATAoffline) ;
%     elseif DATAoffline.ReactiveModesFromInterfaceDisplacements == 1
%         % TENTATIVE VERSION 1-APRIL-2023
%         % WE FOUND IN
%         % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/06_3D_27_HEXAHEDRA/04_FinalExam.mlx
%         % that the version used in giuliudory2023multiscale.pdf is prone to
%         % instabilities. The version implemented here attempts to eliminate
%         % such instabilities
%         [BasisUdeform,PsiDEFf,Vall,Mintf,Mintf_INV,Mdom,PhiRB,PsiRBf,MESHdom,DATA]  =...
%             MODESandROMoper_EIFE(OPERFE,MESHdom,DATA,SNAPreact,BasisStwo,...
%             SNAPdisp,DATAcommon,DATAoffline) ;
%
%     else
%         error('Option not implemented')
%     end
% end



% 5) HYPERREDUCTION (SELECTION POINTS AND WEIGHTS)
%-------------------

%----------------------------------
DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',[]) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'INTERNAL_FORCES',0) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'MASS_MATRIX',0) ;
DATA.BasisSTRESS_SINGULAR_VALUES = BasisSTRESS_SINGULAR_VALUES ;

CECM= ContinuousECM_multi(OPERFE,[PhiDEF,GammaBUB],BasisSTRESS,DATA,DATAcommon,MESHdom,DATAoffline,MATPRO,PhiRB) ;

MODES.PhiDEF = PhiDEF;
MODES.GammaBUB = GammaBUB;
MODES.PsiDEFf = PsiSEf;
MODES.Vall = Vall;
MODES.PhiRB = PhiRB;
MODES.PsiRBf = PsiRBf;
MODES.BasisPone = BasisSTRESS;
MODES.BasisStwo = BasisSTRESS;


GEOMATRICES.Mintf = Mintf;
%GEOMATRICES.MdomFF = MdomFF;
GEOMATRICES.Mdom = Mdom;

DATA_misc = DATA;



save(NAMEOFFLINEstore,'MODES','GEOMATRICES','MESHdom','DATA_misc','DATAoffline','MATPRO','OPERFE','CECM','DATAcommon')



diary off


