function OFFLINE_STEPSbubCOMPL(TrainingTrajectoryINP,DATAoffline)
%--------------------------------------------------------------------------
% FUNCTION: OFFLINE_STEPSbubCOMPL
%
% PURPOSE:
%   Executes the complete offline stage of the multiscale EIFEM methodology
%   using bubble-compliant modes. It constructs the reduced-order basis
%   for displacements and stresses from training data, and applies
%   hyperreduction (Discrete Empirical Cubature Method, CECM) to compress
%   internal force computations.
%
%   Key steps include:
%     - Loading displacement snapshots and mesh/domain data
%     - Constructing the global deformation modes (excluding rigid-body parts)
%     - Projecting stresses and constructing stress basis
%     - Generating bubble modes and self-equilibrated fluxes
%     - Assembling final reduced spaces and applying CECM
%
% INPUTS:
%   - TrainingTrajectoryINP : String name of a function defining
%                             the training configuration (mesh, BCs, etc.).
%
%   - DATAoffline           : (Optional) Struct with parameters for:
%                               * SVD tolerances
%                               * Domain/faces of interest
%                               * Stress mode selection method
%                               * Use of DECM (Discrete ECM) for hyperreduction
%                               * Accuracy checks and storage flags
%
% OUTPUT:
%   - A binary .mat file named [SNAPSHOTS/NAME_BASE_OFFLINE.mat] containing:
%       * TRICES (mass matrices on domain/interface)
%       * OPERFEMODES (RB, deformation, bubble, stress)
%       * GEOMA (FE operators)
%       * MATPRO (material model data)
%       * MESHdom (mesh of the domain of interest)
%       * CECM (cubature weights and reduced points)
%       * DATA_misc, DATAoffline, DATAcommon (metadata and configuration)
%
% RELATED FUNCTIONS:
%   - `GetDisplacAndMesh_sevDOM` : Retrieves mesh and training displacements
%   - `DeformModesNONL_INVAR`    : Constructs invariant deformation modes
%   - `RetrievingSTRESSESprojected` : Projects PK1 stresses to reduced space
%   - `BubbleModes_BasicSEfromDEF` : Computes compliant bubble and SE fluxes
%   - `StressModes_viaSVD`, `StressModes_viaSVDinelCOMPL` : Stress basis
%   - `ContinuousECM_multi` : Computes reduced cubature (CECM)
%
% SEE ALSO:
%   - /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/
%     TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/05_Assessment.mlx
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   Created: 30-Jan-2023 — Last Updated: 23-Sep-2023
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

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

%
% DATAOUT.LABEL_DOMAIN = 1;
% DATAOUT.LABELS_FACES =   [3,4,5,6,7,8] ;   % LABELS IDENTIFYING THE BOUNDARY OF  the  DOMAIN, in the PRE-GID FILE
%                                    % IT INCLUDES INTERFACE BOUNDARIES AND
%                                    % NON-INTERFACE BOUNDARIES
DATAcommon = DefaultField(DATAcommon,'LABEL_DOMAIN',1) ;   % 23-Feb-2025
DATAcommon = DefaultField(DATAcommon,'LABELS_FACES',[]) ;
DATAoffline = DefaultField(DATAoffline,'LABEL_DOMAIN',DATAcommon.LABEL_DOMAIN) ;
DATAoffline = DefaultField(DATAoffline,'LABELS_FACES',DATAcommon.LABELS_FACES) ;

if isempty(DATAoffline.LABELS_FACES)
    disp('You must specify whicha are the boundaries of the selected domain')
end


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


% -------------------------------
% 3) Interface displacements
% --------------------------------

DATAcommon = DefaultField(DATAcommon,'ORIGIN_COOR_CENTROID_INTERFACE',0) ;

if  DATAcommon.ORIGIN_COOR_CENTROID_INTERFACE == 1
    % Local coordinates coincides with the centroid of the interface
    % (Feb-2025),
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
    % boundaries
    [PhiRB,PsiRBf,Mdom,Mintf,MESHdom] =   GeometricVarDOMAINScINTF(OPERFE,MESHdom,DATA,DATAcommon);
else
    % Local coordinates coincides with the centroid of  the domain
    [faceDOFS,PhiRB,PsiRBf,Mdom,Mintf,MESHdom] =   GeometricVarDOMAINS(OPERFE,MESHdom,DATA,DATAcommon);
end

DATA.NAME_BASE = NAME_BASE;
[Vall,MESHdom,DATA] = FictInterfaceDISP_raw(MESHdom,DATAcommon,DATA,Mintf,DATAoffline) ;



% 4) Compute LINEAR and NONLINEAR deformation modes (basic + complementary) 
% using SVD and mass-orthogonalization; also returns stiffness and GEOMETRIC mass MATRICES


[BasisUdeform,MESHdom,DATA,DATAoffline,Kstiff,MdomCHOL,MintfINV_chol]  =...
    DeformModesNONL_INVAR(OPERFE,MESHdom,DATA,...
    SNAPdisp,DATAcommon,DATAoffline,Vall,PhiRB,Mdom,Mintf,MATPRO,SNAPdisp_AUX) ;

 

% -------------------------------------------
% 6) RETRIEVING STRESSES AND REACTIVE FORCES  (project-wise)
% ------------------------------------------
%DATAoffline = DefaultField(DATAoffline,'CECM_ONLY_FOR_NONLINEAR_STRESSES',1) ; % 6-Jan-2023,  see
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
[AsnapSTRESS_PK1,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel] ...
    = RetrievingSTRESSESprojected(DATAoffline,NAME_BASE,DATAcommon,NAMEsnap_base,INFO_RVE,BasisUdeform,DATA,...
    OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;

% ------------------------------------------------
% 3) Self-equilibrated modes and Bubble Modes
% ------------------------------------------------
DUMMY = [] ;
[PsiSEf,PhiDEF,GammaBUB] =     BubbleModes_BasicSEfromDEF(PsiRBf,Mintf,DUMMY,DATA,DATAoffline,MESHdom,...
    BasisUdeform,Mdom,Kstiff,MATPRO,OPERFE,MintfINV_chol,PhiRB,MdomCHOL)  ;

% ----------------------------------------------------------------------------------------------
% 4) Check accuracy after bubble modes
BasisUdeform_bubble = [PhiDEF,GammaBUB] ;
CheckAccuracyStressesAfterBUB(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
    INFO_RVE,BasisUdeform_bubble,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB)
% ------------------------------------------------------------------------------------------------

% 4) Stress modes
% ---------------------------------------------------------------------


if ~isempty(AsnapSTRESSel)
    [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_viaSVDinelCOMPL(MATPRO,AsnapSTRESSinel,DATA,DATAoffline,MESHdom,OPERFE,BasisUdeform,AsnapSTRESSel)  ;
else
    [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] =   StressModes_viaSVD(AsnapSTRESS_PK1,DATA,DATAoffline,MESHdom,OPERFE,BasisUdeform)  ;
    
end


%


% CASE OF 1D ELEMENTS
% --------------------
switch DATAcommon.TypeFunctionDisplacementInterfaces
    case '1D_LINEAR'
        disp('1D elements...')
        disp('Only one rigid-body mode...')
        PhiRB = PhiRB(:,1) ;
        PsiRBf = PsiRBf(:,1) ;
        
end


% 5) HYPERREDUCTION (SELECTION POINTS AND WEIGHTS)
%-------------------

%----------------------------------
DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',[]) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'INTERNAL_FORCES',0) ;
DATAoffline.UseDECMpoints = DefaultField(DATAoffline.UseDECMpoints,'MASS_MATRIX',0) ;
DATA.BasisSTRESS_SINGULAR_VALUES = BasisSTRESS_SINGULAR_VALUES ;

DATA =DefaultField(DATA,'INDEXES_LINEAR_MODES_DEFORMATION',[]) ;
%
% CONSIDER_SEPARATED_SVDsINTERNALFORCES = 0 ;
%
% if
%
% if ~isempty(DATA.INDEXES_LINEAR_MODES_DEFORMATION)
%     if ~isempty(DATA.INDEXES_LINEAR_MODES_DEFORMATION.COMPLEMENTARY)
%         PhiDEFall = cell(1,2) ;
%         PhiDEFall{1}
%     else
%         PhiDEFall = [PhiDEF,GammaBUB] ;
%     end
%
% else
%
%
% end
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



save(NAMEOFFLINEstore,'MODES','GEOMATRICES','MESHdom','DATA_misc','DATAoffline','MATPRO','OPERFE','CECM','DATAcommon','Kstiff')



diary off


