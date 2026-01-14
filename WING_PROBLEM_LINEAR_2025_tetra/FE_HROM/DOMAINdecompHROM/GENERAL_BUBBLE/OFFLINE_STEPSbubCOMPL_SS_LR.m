% function OFFLINE_STEPSbubCOMPL_SS_LR(TrainingTrajectoryINP,DATAoffline)
% % Modified version of OFFLINE_STEPSbubCOMPL.m, to be invoked from 
% %  Train1dom_BUBcomp_SS_LR.m
% % It covers scenarions in which either strains are small, but rotations
% % large, or strains are large, but rotations small. 
% 
% % Offline steps for constructing a hyperreduced-order  MULTISCALE MODEL
% % Bubble modes
% % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
% % JAHO, previous versions 30-JAN-2023/23-Sept-2023/17-Oct-2023
% % Created on 4-Feb-2025, Balmes 185, Barcelona. 
% %-------------------------------------------------------------
% if nargin == 0
%     load('tmp.mat')
% end
% delete('OFFLINE.txt')
% diary 'OFFLINE.txt'
% if ~exist('Quadratic3N')
%     addpath(genpath(getenv('MATLAB_CODES_FEHROM')))  ; %
% end
% % for controling SVD accuracy
% % DATAoffline = DefaultField(DATAoffline,'errorDISP',0) ;
% % DATAoffline = DefaultField(DATAoffline,'errorSTRESS',0) ;
% DATAoffline = DefaultField(DATAoffline,'errorFINT',1e-10) ;
% % DATAoffline = DefaultField(DATAoffline,'errorECM',0) ;
% % DATAoffline = DefaultField(DATAoffline,'errorPK2stress_basis',DATAoffline.errorSTRESS) ;
% % DATAoffline = DefaultField(DATAoffline,'errorREACTIVE_FORCES',0) ;
% % DATAoffline = DefaultField(DATAoffline,'errorSVD_SELF_EQUILIBRATED',1e-6) ;
% 
% DATAoffline = DefaultField(DATAoffline,'TOLSVD_complementary_modes',1e-3) ;
% 
% 
% 
% DATAoffline = DefaultField(DATAoffline,'TOLSVD_complementary_modes_stresses',1e-10) ;
% DATAoffline = DefaultField(DATAoffline,'BASIS_STRESSES_FROM_PROJECTED_DISPLACEMENTS',1) ;
% 
% 
% % DATAoffline.TOLSVD_complementary_modes = 1e-3 ; % displacements
% % DATAoffline.TOLSVD_complementary_modes_stresses = 1e-3 ;
% % DATAoffline.BASIS_STRESSES_FROM_PROJECTED_DISPLACEMENTS = 0 ;
% 
% 
% 
% %
% % READING INPUT DATA
% DATAcommon = feval(TrainingTrajectoryINP) ;
% DATAcommon = DefaultField(DATAcommon,'NameFileMeshDATA',[] ) ;
% DATAoffline = DefaultField(DATAoffline,'NAMEMESH', DATAcommon.NameFileMeshDATA);
% DATAcommon.NameFileMeshDATA = DATAoffline.NAMEMESH ;
% 
% NAME_BASE =[DATAcommon.NameParamStudyLOC,'_param_'];
% NAMEsnap_base = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE];
% NAMEOFFLINEstore = [cd,filesep,'SNAPSHOTS',filesep,NAME_BASE,'OFFLINE.mat'];
% 
% % -------------------------------------------
% % 1) RETRIEVING DISPLACEMENT MODES AND MESH INFORMATION
% % ------------------------------------------
% DATAoffline.ReturnFullDisplacementMatrix = 1;
% DATAoffline.USE_BASIS_MATRIX_DISPLACEMENTS_AS_SNAPSHOTS = 2;
% [INFO_RVE,MESHdom,SNAPdisp,DATA,OTHER_output,OPERFEaux,MATPRO,Fbody,Ftrac,SNAPdisp_AUX,DATAoffline] ...
%     = GetDisplacAndMesh_sevDOM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline) ;
% 
% 
% % -------------------------------------------------------------------------------------------------------------------------------------------------
% % 2) FE OPERATORS AND OTHER VARIABLES RELATED WITH THE FE DISCRETIZATION OF THE DOMAIN (SUCH AS GLOBAL SHAPE FUNCTIONS, Global B-matrix operator)
% % --------------------------------------------------------------------------------------------------------------------------------------------------
% [OPERFE,MESHdom,DATA] = FEoperators_1dom(DATA,MESHdom,MATPRO)   ;
% 
% 
% % -------------------------------
% % 3) Interface displacements
% % --------------------------------
% 
% DATAoffline.ORIGIN_COOR_CENTROID_INTERFACE = 1; % THIS IS THE DEFAULT OPTION ! 
% 
% %if  DATAoffline.ORIGIN_COOR_CENTROID_INTERFACE == 1
%      [PhiRB,PsiRBf,Mdom,Mintf,MESHdom] =   GeometricVarDOMAINScINTF(OPERFE,MESHdom,DATA,DATAcommon);
% % else
% %     [faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESHdom] =   GeometricVarDOMAINS(OPERFE,MESHdom,DATA,DATAcommon);
% % end
% 
% DATA.NAME_BASE = NAME_BASE;
% [Vall,MESHdom,DATA] = FictInterfaceDISP_raw(MESHdom,DATAcommon,DATA,Mintf,DATAoffline) ;
% 
% 
% 
% 
% % DEFORMATIONAL MODES (PURGE RIGID-BODY COMPONENT).
% % -------------------
% 
% %DATAoffline = DefaultField(DATAoffline,'METHOD_INVARIANT_DEF_MODES',[]);
% 
% %DATAoffline.METHOD_INVARIANT_DEF_MODES = DefaultField(DATAoffline.METHOD_INVARIANT_DEF_MODES,'STAGE_TRAINING',1);
% 
% DATAcommon = DefaultField(DATAcommon,'METHOD_TRAINING_ELASTIC_RANGE','INVARIANT_MODES') ;
% 
% 
% if  ~isempty(DATAcommon.METHOD_TRAINING_ELASTIC_RANGE)
%     % ------------------------------------------------
%     % Specific method for determining invariant modes + complementary
%     % modes,
%     % Apr-3-2024
%     % see /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/107_EIFEM_plast3D/02_BEAMQ8_IprofE.mlx
%     % -----------------------------------------------
%     switch  DATAcommon.METHOD_TRAINING_ELASTIC_RANGE
%         case  'INVARIANT_MODES_PLUS_COMPLEMENTARY_MODES'
%             disp(' Not used....Try better "invariant_modes"')
%             [BasisUdeform,MESHdom,DATA,DATAoffline,Kstiff]  =...
%                 DeformModesNONL_INVARIANTSC(OPERFE,MESHdom,DATA,...
%                 SNAPdisp,DATAcommon,DATAoffline,Vall,PhiRB,Mdom,Mintf,MATPRO) ;
%             MdomCHOL=[] ; MintfINV_chol = [] ;
%             
%         case 'INVARIANT_MODES'
%             [BasisUdeform,MESHdom,DATA,DATAoffline,Kstiff,MdomCHOL,MintfINV_chol]  =...
%                 DeformModesNONL_INVAR(OPERFE,MESHdom,DATA,...
%                 SNAPdisp,DATAcommon,DATAoffline,Vall,PhiRB,Mdom,Mintf,MATPRO) ;
%             
%             
%         otherwise
%             error('Method not implemented')
%             
%     end
%     
% else
%     error('Option not implemented in this OFFLINE_STEP file')
%     % See rather OFFLINE_STEPbub.m
% end
% 
% 
% % -------------------------------------------
% % 2) RETRIEVING STRESSES AND REACTIVE FORCES  (project-wise)
% % ------------------------------------------
% %DATAoffline = DefaultField(DATAoffline,'CECM_ONLY_FOR_NONLINEAR_STRESSES',1) ; % 6-Jan-2023,  see
% % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
% 
% 
% 
% 
% DATAoffline = DefaultField(DATAoffline,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; % 6-Jan-2023,  see
% % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
% 
% if DATAoffline.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
%     AsnapSTRESSinel = [] ;
%     AsnapSTRESSel =[] ;
%     [AsnapREAC ,AsnapSTRESS_PK2 ,AsnapSTRESS_PK1,NAME_BASE]  = GetStressesAndReactForces_bub(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
%         INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
% else
%     % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
%     % Separate treatment nonlinear stresses
%     % Separate treatment nonlinear stresses
%     
%     if  DATA.SMALL_STRAIN_KINEMATICS ==1  && DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1
%         [AsnapREAC ,AsnapSTRESS_PK2 ,AsnapSTRESS_PK1,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel]  =...
%             GetStressesAndReactForces_bubNECM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
%             INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
%         
%     else
%         [AsnapREAC ,AsnapSTRESS_PK2 ,AsnapSTRESS_PK1,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel] ...
%             = GetStressesAndReactForces_bubNECMlarg(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
%             INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
%         
%     end
%     
% end
% 
% 
% 
% 
% 
% % 3) Self-equilibrated modes and Bubble Modes
% %[PsiSEf,PhiDEF,GammaBUB] =     BubbleModes_BasicSE(PsiRBf,Mintf,AsnapREAC,DATA,DATAoffline,MESHdom,BasisUdeform,Mdom)  ;
% 
% [PsiSEf,PhiDEF,GammaBUB] =     BubbleModes_BasicSEfromDEF(PsiRBf,Mintf,AsnapREAC,DATA,DATAoffline,MESHdom,...
%     BasisUdeform,Mdom,Kstiff,MATPRO,OPERFE,MintfINV_chol)  ;
% 
% DATAoffline =DefaultField(DATAoffline,'CHECK_ACCURACY_PROJECTED_STRESSES_AFTER_BUBBLE_MODES',0) ;
% 
% if DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES_AFTER_BUBBLE_MODES == 1
%     
%     BasisUdeform_bubble.PhiDEFbs = PhiDEF ;
%     BasisUdeform_bubble.PhiDEFcomp = GammaBUB;
%     
%     
%     DATAoffline.CHECK_ACCURACY_PROJECTED_STRESSES = 1 ;
%     [~ ,~ ,~,~]  = GetStressesAndReactForces_bub(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
%         INFO_RVE,BasisUdeform_bubble,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
%     
% end
% 
% 
% 
% % 4) Stress modes
% % ---------------------------------------------------------------------
% 
% 
% if ~isempty(AsnapSTRESSel)
%     
%     [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_viaSVDinelCOMPL(MATPRO,AsnapSTRESSinel,DATA,DATAoffline,MESHdom,OPERFE,BasisUdeform,AsnapSTRESSel)  ;
%     
% else
%     [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_viaSVD(AsnapSTRESS_PK1,DATA,DATAoffline,MESHdom,OPERFE,BasisUdeform)  ;
%     
% end
% 
% 
% %
% 
% 
% % CASE OF 1D ELEMENTS
% % --------------------
% switch DATAcommon.TypeFunctionDisplacementInterfaces
%     case '1D_LINEAR'
%         disp('1D elements...')
%         disp('Only one rigid-body mode...')
%         PhiRB = PhiRB(:,1) ;
%         PsiRBf = PsiRBf(:,1) ;
%         
% end
% 
% 
% % 5) HYPERREDUCTION (SELECTION POINTS AND WEIGHTS)
% %-------------------
% 
% %----------------------------------
% DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',[]) ;
% DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'INTERNAL_FORCES',0) ;
% DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'MASS_MATRIX',0) ;
% DATA.BasisSTRESS_SINGULAR_VALUES = BasisSTRESS_SINGULAR_VALUES ;
% 
% 
% CECM= ContinuousECM_multi(OPERFE,[PhiDEF,GammaBUB],BasisSTRESS,DATA,DATAcommon,MESHdom,DATAoffline,MATPRO,PhiRB) ;
% 
% MODES.PhiDEF = PhiDEF;
% MODES.GammaBUB = GammaBUB;
% MODES.PsiDEFf = PsiSEf;
% MODES.Vall = Vall;
% MODES.PhiRB = PhiRB;
% MODES.PsiRBf = PsiRBf;
% MODES.BasisPone = BasisSTRESS;
% MODES.BasisStwo = BasisSTRESS;
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
% save(NAMEOFFLINEstore,'MODES','GEOMATRICES','MESHdom','DATA_misc','DATAoffline','MATPRO','OPERFE','CECM','DATAcommon','Kstiff')
% 
% 
% 
% diary off
% 
% 
