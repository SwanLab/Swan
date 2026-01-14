function  [BasisSTRESS,BasisSTRESS_SINGULAR_VALUES ] = ...
    StressModes_viaSVDinelCOMPL(MATPRO,AsnapSTRESS,DATALOC,DATAoffline,MESH,OPERFE,BasisUdeform,AsnapSTRESSel)  ;

%==========================================================================
% function [BasisSTRESS, BasisSTRESS_SINGULAR_VALUES] = ...
%   StressModes_viaSVDinelCOMPL(MATPRO, AsnapSTRESS, DATALOC, ...
%                                DATAoffline, MESH, OPERFE, ...
%                                BasisUdeform, AsnapSTRESSel)
%
% PURPOSE:
% Extracts a reduced basis of **inelastic stress modes** via a weighted
% SVD (WSVD) of stress snapshots collected during nonlinear simulations.
% These stress modes serve as the second-layer reduction basis for the
% EIFEM method (Empirical Interscale FEM), particularly in the construction
% of CECM-based hyperreduced models.
%
% The focus is on the **complementary** (i.e., nonlinear or inelastic)
% stress contributions, discarding basic (elastic) components unless
% explicitly required (e.g., when no nonlinear stresses are available).
%
%==========================================================================
% THEORETICAL CONTEXT:
% The procedure implemented follows the theoretical description in
% Sections "Hyperreduction of nodal internal forces" and
% "Efficient computation of the snapshot matrix of the integrand"
% (pp. 222–235, ). The approach can be summarized as follows:
%
% 1. **Stress Snapshots Collection:**
%    - From nonlinear training tests, collect the stress field at all Gauss
%      points for each time step: ∆Pf′(Xg ; ti)
%    - Store these as columns in a matrix SΔP
%
% 2. **Weighting:**
%    - A diagonal weighting matrix W is defined using the interpolation
%      weights (typically the product of Gauss weights and Jacobians)
%    - This ensures SVD extracts modes with respect to the L2-inner product
%
% 3. **Compression via WSVD:**
%    - A (possibly truncated) weighted SVD is applied:
%        [Λ, Σ, V] ← WSVD(SΔP, W)
%    - The resulting left singular vectors Λ define the inelastic stress modes
%
% 4. **Fallback Mechanism:**
%    - If inelastic snapshots are unavailable, elastic stress snapshots
%      (PK1 or small-strain) are synthesized using:
%         σ = C : ε(BΦ)
%      where B is the strain-displacement matrix and Φ are the deformational modes
%
%==========================================================================
% INPUTS:
% - MATPRO          : Structure with material stiffness matrix (Celas or celasglo)
% - AsnapSTRESS     : Cell array of stress snapshots (both linear and nonlinear)
% - DATALOC         : Local configuration parameters (e.g., strain type)
% - DATAoffline     : Structure containing tolerances and settings
% - MESH            : Mesh structure (provides #strains, DOFs, etc.)
% - OPERFE          : Contains operator wSTs for weighting
% - BasisUdeform    : Deformational basis used to construct stresses
% - AsnapSTRESSel   : Elastic stress snapshots (optional fallback)
%
% OUTPUTS:
% - BasisSTRESS              : Matrix of stress modes extracted via WSVD
% - BasisSTRESS_SINGULAR_VALUES : Singular values corresponding to each mode
%
%==========================================================================
% NOTES:
% - This function is typically called after bubble modes have been computed.
% - The stress modes are used to construct reduced-order internal force
%   models for hyperreduction (CECM).
% - When nonlinear tests are unavailable, a warning is displayed and
%   elastic approximations are used (e.g., PK1 stress from Φ).
%
%==========================================================================
% AUTHOR:
% Joaquín Alberto Hernández Ortega, Campus Nord UPC, Barcelona
% DATE:
% Original: 08-Apr-2024
% Annotated and commented: May 2025 (ChatGPT)
%==========================================================================


% ------------------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
    %  DATAoffline.TOLSVD_complementary_modes = 1e-3 ;
end


% Basic stress modes
% -----------------

% THERE SHOULD BE NO BASIC STRESS MODES HERE
% disp('-----------------------------------------')
% disp('Basis STRESS modes')
% disp('-----------------------------------------')
%AsnapSTRESS_basic = cell2mat(AsnapSTRESS(DATAoffline.INDEX_BASIC_TESTS)) ;
%
%  DATALOC.LEGEND_MODES_stress = 'OriginalSnapshotsBasic' ;
%
% PlotModesStress_SUBDOMAIN(DATALOC,AsnapSTRESS_basic,MESH);


% Basis matrix for
% DATALOC = DefaultField(DATALOC,'TOL',1e-6);
% % -------------------------------------
%
% % wST = repmat()
%
% % A = spdiags(Bin,d,6,6);
%
% % wST = speye(length(OPERFE.wSTs) ; % Diagonal matrix with weights

% wST = speye(length(OPERFE.wSTs) ; % Diagonal matrix with weights

%  COMPLEMENTARY MODES STRESSES
disp('-----------------------------------------')
disp('STRESS MODES INELASTIC STRESSES' )
disp('-----------------------------------------')

AsnapSTRESS_compl= cell2mat(AsnapSTRESS(DATAoffline.INDEX_COMPL_TESTS)) ;  % THESE ARE

nstrain = DATALOC.MESH.nstrain ;
ISLARGE = 0 ;
if  DATALOC.NO_USE_Deformation_gradient_in_Small_Strains ==0 || DATALOC.SMALL_STRAIN_KINEMATICS == 0
    % Large strain formulation
    nstrain = DATALOC.MESH.ndim^2;
    ISLARGE = 1;
end

wST = repmat(OPERFE.wSTs',nstrain,1);
wST = wST(:) ;
W = spdiags(wST,0,length(wST),length(wST));
%
%




if isempty(AsnapSTRESS_compl)
    disp(['Problem with no inelastic stresses'])
    disp(['Using elastic stresses for computing CECM points (use low tolerance for internal forces if you want to reduce the number of points)'])
    
    if ISLARGE == 1
        % JAHO, 24-MAY-2024
        % CASE LARGE STRAINS. ELASTIC STRESSES (PK1) ARE COMPUTED IN 
        % MultiSnapStressFromDispNECMlarg.m
        % See 
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
        AsnapSTRESS_compl  = cell2mat(AsnapSTRESSel) ; 
         
    else
        
        if isfield(MATPRO,'Celas')
            MATPRO.celasglo = MATPRO.Celas ;
        end
        
        
        if ~isempty(MATPRO.celasglo)
            %     disp('Checking accuracy CECM points ')
            %     disp('-------------------------------------------')
            %   disp('Coarse-scale stiffness matrix ')
            %     celastST = MATPRO.celasglo ;
            %     nF = size(MATPRO.celasglo,2) ;
            %     for icomp = 1:nF
            %         icol = icomp:nF:size(celastST,1) ;
            %         celastST(icol,:) = bsxfun(@times,celastST(icol,:),OPERFE.wSTs) ;
            %     end
            celastST = ConvertBlockDiag(MATPRO.celasglo) ; % Diagonal block matrix
            AsnapSTRESS_compl =  celastST*(OPERFE.Bst*BasisUdeform.PhiDEFbs);
        else
            error('This option requires MATPRO.celasglo ')
        end
        
    end
    
    
    
end



DATALOC = DefaultField(DATALOC,'TOL',1e-6);



Wchol = sqrt(W);
A = {Wchol*AsnapSTRESS_compl} ;
TOL = DATAoffline.TOLSVD_complementary_modes_stresses; % Tolerance second block
RELTOL = [TOL] ;
DATAlocS = [] ;
[U,S,V] = SRSVD(A,RELTOL,DATAlocS) ;

U = Wchol\U ;  % All modes (including basic ones)

BasisSTRESS = U ;
BasisSTRESS_SINGULAR_VALUES = S ;
disp(['Number of nonlinear stress modes =',num2str(length(S))])

disp('Singular values')

S

 
disp('PLotting stress modes...(nonlinear contribution)')
if ISLARGE == 1
    disp('Warning = PK1 stresses are not correctly post-process in GID (they are not symmetric)')
    disp('Se we cannot plot the nonlinear contribution of the stresses')
    
else
    DATALOC.LEGEND_MODES_stress = 'Inelastic_stress ()' ;
    
    PlotModesStress_SUBDOMAIN(DATALOC,BasisSTRESS,MESH);
    
end





% JAHO 22-Apr-2024, to avoid large matrices of internal forces, if not
% required
DATAoffline = DefaultField(DATAoffline,'NumberOfStressModesGiven',[]) ;
if ~isempty(DATAoffline.NumberOfStressModesGiven)
    NumberOfStressModesGiven = min(DATAoffline.NumberOfStressModesGiven,length(S)) ;
    
    BasisSTRESS = U(:,1:NumberOfStressModesGiven) ;
    BasisSTRESS_SINGULAR_VALUES = S(NumberOfStressModesGiven) ;
end



% -------------------------------------------------
% POST-PROCESSING SELF-EQUILIBRATED MODES
% -------------------------------------------------
