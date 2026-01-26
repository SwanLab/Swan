function  [xNEW,wNEW,VAR_SMOOTH_FE,DATAIN,DATA_REFMESH] = OptimalECMpoints(DATAIN,DATA_GENGAUSS,DATA_REFMESH,DATAROM,BASES)
% This algorithm returns a set of optimal integration points
% It is an extension to continuous  FE structures of the Empirical Cubature
% Method.
% INPUTS
% ------
% DATAIN --> Sundry FE inputs (loaded from FE results)
% DATA_GENGAUSS --> Input data specific to the optimization problem
% DATA_REFMESH --> FE inputs (COORDINATES, connectivities...)
% BASES, DATAROM ---> Inputs related to the FE reduced-order model (basis matrices, etc...)

% OUTPUTS
% --------
% xNEW --> coordinates optimal points
% wNEW --> Weights optimal points
% VAR_SMOOTH_FE --> Sundry variables (some of them are  outputs of the smoothing process)
% DATAIN ...

% April 26th 2020. JAHO (44th day of confinement, COVID 19)

DATAOUT = [] ;
MSG = {} ;

switch   DATA_GENGAUSS.APPROX_FUN__DERI.METHOD
    case 'FE_ONLY_NOdecomposition'
        % Disabling some options (no stress and B matrices in this method)
        DATA_GENGAUSS.PLOT_Basis_STRESS  = 0 ;
        DATA_GENGAUSS.PLOT_Basis_Bred_MODES  = 0 ;
        error('Option not implemented yet')
end


% Checking incomp.
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'TOL_LOC_InternalForces_IS_INTEGRATION',0) ;  % 9th-April-2020 (27th day quarantine COVID19)
DATAIN.TOL_LOC_InternalForces_IS_INTEGRATION = DATA_GENGAUSS.TOL_LOC_InternalForces_IS_INTEGRATION; % (Not relevant 26-4-20)
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'APPROX_FUN__DERI',[]) ;  %
DATA_GENGAUSS.APPROX_FUN__DERI=DefaultField(DATA_GENGAUSS.APPROX_FUN__DERI,'METHOD','FE_ONLY'); % 26-4-20 % Default Value from 26-4-20
DATAIN  = DefaultField(DATAIN,'CUBATURE',[]) ;
DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'IMPOSE_VOLUME_CONSTRAINT',0) ; % Default value = 0. In FE multiscale problems, the volume is exactly integrated
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'TOL_LOC_InternalForces_SVD',1e-6) ;  %  Tolerance SVD internal forces matrix
DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'TOL_ECM',0) ;  %  Toleranc ECM

DATAIN.CUBATURE.TOL_LOC_InternalForces  = DATA_GENGAUSS.TOL_LOC_InternalForces_SVD ;

% Selecting candidates for the ECM
% ----------------------------------------------------------------
% Here  we have the possibility of selecting which points of the
% underlying FE mesh we  wish to include in the set of candidate ECM points
%
% -----------------------------------------------------------------------------------------------------------
DATA_REFMESH = DefaultField(DATA_REFMESH,'DOFr',[]) ; % Slave/Constrained DOFs
DATA_GENGAUSS.TypeElement = DATA_REFMESH.TypeElement  ; 

[INDSEL,DATAdistorsion ]= RestrictedDomainForECMpoints(DATA_GENGAUSS,DATA_REFMESH.COOR,DATA_REFMESH.Nst,DATA_REFMESH.CN,...
    DATA_REFMESH.CONNECTb) ;
MSG{end+1}  = DATAdistorsion.MSG ; 
DATA_REFMESH.DATAdistorsionELEM = DATAdistorsion; 
%------------------------------------------------------------------------------------------------------------
DATAIN.GAUSS_POINTS_TO_BE_INCLUCED_IN_THE_ECM = INDSEL ;
DATAIN.CUBATURE.TOLERANCE_EmpiricalCubature =  DATA_GENGAUSS.TOL_ECM;

% Deteermine ECM points
[HROMVAR,MSG,wSTs,BasisS_gauss,BdomRED_gauss,SingVal_F,VrightVal_F,DATA_REFMESH,SNAPforceSnw]= ...
    ReducedSetIntegrationPoints_loc(DATAIN,DATAROM,BASES,DATA_REFMESH,MSG,DATA_GENGAUSS) ;


% ------------------------------------------------------------------------------------------
switch DATA_GENGAUSS.APPROX_FUN__DERI.METHOD
    case {'FE_INTERPOLATION','FE_ONLY'}
        % Computing stress and virtual strains basis matrices at the nodes
        % of the FE discretization
        [VAR_SMOOTH_FE,MSG] = ...
            BasisAtNodes_fromGAUSS(DATA_GENGAUSS,DATA_GENGAUSS.NAME_INPUT_DATA,...
            BasisS_gauss,BdomRED_gauss,DATAIN,DATA_REFMESH,HROMVAR,MSG,SingVal_F,VrightVal_F,wSTs) ;
        % Variables required for computing the Jacobian Matrix in the
        % removal iterative proceudre
        
    otherwise
        BdomRED_nodes = [] ; BasisS_nodes = [] ;
        
end
VAR_SMOOTH_FE.wSTs =  wSTs ;
VAR_SMOOTH_FE.IMPOSE_VOLUME_CONSTRAINT = DATA_GENGAUSS.IMPOSE_VOLUME_CONSTRAINT ;

DATAIN.nstrain = HROMVAR.nstrain ;



%%% POINTS REQUIRED TO INTEGRATE STIFFNESS MATRIX
% -----------------------------------------------
% disp('-----------------------------------------------------')
% disp('POINTS REQUIRED TO INTEGRATE STIFFNESS MATRIX')
% disp('-----------------------------------------------------')
DATAIN = DefaultField(DATAIN,'MSGPRINT',{}) ;
for iii = 1:length(MSG)
    DATAIN.MSGPRINT{end+1} = MSG{iii} ;  
end
 
% Checking integration error
% -----------------------------
EXACTint = HROMVAR.PHI'*wSTs ;
APPROXint = HROMVAR.PHI(HROMVAR.setPoints,:)'*HROMVAR.WdomRED ;
eee = norm(EXACTint-APPROXint)/norm(EXACTint)*100 ;
DATAIN.MSGPRINT{end+1} = ['Difference between approx. and ex. int  of LEFT SINGULAR VECTORS =',num2str(eee),' %'] ;
disp(DATAIN.MSGPRINT{end} )


 

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_INTERNAL_FORCE_MODES',0) ;

if DATA_GENGAUSS.PLOT_INTERNAL_FORCE_MODES == 1
    DATAIN =   GidPostProcessModesFINT_loc(DATA_REFMESH.COOR,DATA_REFMESH.CN,...
        DATA_REFMESH.TypeElement,HROMVAR.PHI,DATA_REFMESH.posgp,DATA_GENGAUSS.NAME_INPUT_DATA,DATAIN);
end

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_INTERNAL_FORCE_SNAPSHOTS',0) ;  %

if DATA_GENGAUSS.PLOT_INTERNAL_FORCE_SNAPSHOTS == 1
SNAPforceSnw_N  =norm(SNAPforceSnw,'fro') ; 
SNAPforceSnw = SNAPforceSnw/SNAPforceSnw_N ; 
    DATAIN =   GidPostProcessModesFINT_loc(DATA_REFMESH.COOR,DATA_REFMESH.CN,...
        DATA_REFMESH.TypeElement,SNAPforceSnw,DATA_REFMESH.posgp,DATA_GENGAUSS.NAME_INPUT_DATA,DATAIN);
end




DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_Basis_STRESS',0) ;  %
if DATA_GENGAUSS.PLOT_Basis_STRESS == 1  
    NameLOC = [DATA_GENGAUSS.NAME_INPUT_DATA,'_','stressmod'] ;
    MSGLOC = {};
    WHICHMODE = 'STRESS' ;
    DATAstress = [] ;
    MSGLOC = GidPostProcessModes_domGAUSS(DATA_REFMESH.COOR,DATA_REFMESH.CN,...
        DATA_REFMESH.TypeElement,BasisS_gauss,DATA_REFMESH.posgp,NameLOC,DATAstress,WHICHMODE,MSGLOC);
    for iii =1:length(MSGLOC)
        DATAIN.MSGPRINT{end+1} = MSGLOC{iii} ;
    end
end
 

DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_Basis_Bred_MODES',0) ;  %
if DATA_GENGAUSS.PLOT_Basis_Bred_MODES == 1
    NameLOC = [DATA_GENGAUSS.NAME_INPUT_DATA,'_','Bredmod'] ;
    MSGLOC = {};
    WHICHMODE = 'Bred modes' ;
    DATAstress = [] ;
    MSGLOC = GidPostProcessModes_domGAUSS(DATA_REFMESH.COOR,DATA_REFMESH.CN,...
        DATA_REFMESH.TypeElement,BdomRED_gauss,DATA_REFMESH.posgp,NameLOC,DATAstress,WHICHMODE,MSGLOC);
    for iii =1:length(MSGLOC)
        DATAIN.MSGPRINT{end+1} = MSGLOC{iii} ;
    end
end


disp('Generalized Gaussian Cubature. Removing points and adjusting weights')
disp('*************************************************************************++')
DATA_GENGAUSS.posgp = DATA_REFMESH.posgp ;
DATA_REFMESH = DefaultField(DATA_REFMESH,'NameFileMesh',[]) ;
DATA_GENGAUSS.NameFileMesh = DATA_REFMESH.NameFileMesh ;


switch DATA_GENGAUSS.APPROX_FUN__DERI.METHOD
    case {'FE_INTERPOLATION','FE_ONLY'}
        [xNEW,wNEW,DATAIN,ELEMENTS_xNEW,HISTORY_ITERATIONS,VAR_SMOOTH_FE] =  GeneralizedGauss2D_3D_FE...
            (wSTs,DATA_REFMESH.CN,DATA_REFMESH.COOR,DATAIN,...
            DATA_REFMESH.TypeElement,HROMVAR,DATA_REFMESH.Nst,DATA_GENGAUSS,VAR_SMOOTH_FE)  ;
    otherwise
        ELEMENTS_xNEW = [];HISTORY_ITERATIONS = [] ; 
        % All methods implemented before  19th-April-2020
        [xNEW,wNEW,DATAIN] =  GeneralizedGauss2Drectangular...
            (wSTs,DATA_REFMESH.CN,DATA_REFMESH.COOR,DATAIN,...
            DATA_REFMESH.TypeElement,HROMVAR,DATA_REFMESH.Nst,DATA_GENGAUSS,VAR_SMOOTH_FE)  ;
end

VAR_SMOOTH_FE.ELEMENTS_xNEW = ELEMENTS_xNEW ;
VAR_SMOOTH_FE.COORiniECM  = HISTORY_ITERATIONS.POINTS{1}; % iNITIAL SET OF ECM POINTS  
VAR_SMOOTH_FE.WEIGHTSiniECM  = HISTORY_ITERATIONS.WEIGHTS{1}; % iNITIAL SET OF ECM POINTS  
VAR_SMOOTH_FE.ELEMENTSiniECM  = HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{1}; % iNITIAL SET OF ECM POINTS  



disp('*******************************************************************************')

if ~isempty(ELEMENTS_xNEW)
    DATAIN.MSGPRINT{end+1} = 'MESH CONTAINING THE NEW REDUCED SET OF POINTS';
    HROMVAR.WdomRED = wNEW ;
    HROMVAR.setElements = ELEMENTS_xNEW ;
    DATAIN.LABEL_NAME_PROJECT = DATA_GENGAUSS.NAMEFILE ;
    DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_REDUCED_GAUSS_IN_GID',0) ;
    DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'SCALE_FACTOR_FOR_PRINTING_ECM_POINTS',0.1) ;
    
    if DATA_GENGAUSS.PLOT_REDUCED_GAUSS_IN_GID == 1
        DATAIN.PRINT_MESH_ECM_POINTS_AS_POINTS = 1;
        DATAIN.COORDINATES_POINTS_TO_PRINT = xNEW ;
        DATAIN.COORDINATES_POINTS_TO_PRINT = xNEW ;
        DATAIN.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS = DATA_GENGAUSS.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS ;
    end
    
    
    DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS',3) ; % =    1; % 0--> NONE % 1 -->Like time ... % 2 = Separated meshes
    
    if DATA_GENGAUSS.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS == 0
        DATAIN = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HROMVAR) ;
    elseif DATA_GENGAUSS.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS == 1
        DATAIN = PrintingGid_ECMhistory(DATAIN,DATA_REFMESH,HROMVAR,HISTORY_ITERATIONS) ;
    elseif DATA_GENGAUSS.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS >=2
        DATAIN = PrintingGid_ECMpoints_iter(DATAIN,DATA_REFMESH,DATA_GENGAUSS,HISTORY_ITERATIONS) ;
    else
        error('Option not implemented')
    end
    
end