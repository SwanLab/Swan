function [SNAPreactFORCESproj,SNAPstressSTWOproj,SNAPstressPonePROJ,NAME_BASE,SNAPstressPONEproj_NONLIN,SNAPstressPONEproj_LINEAR]  =...
    GetStressesAndReactForces_bubNECMlarg(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
    INFO_RVE,BasisDEFdisp,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB)
% ---------------------------------------------------------------------------------------------------
% FUNCTION: GetStressesAndReactForces_bubNECMlarg
% ---------------------------------------------------------------------------------------------------
% PURPOSE:
%   This function computes **projected first Piola–Kirchhoff (PK1) stresses**, **second Piola–Kirchhoff (PK2) stresses**,
%   and **reaction forces** from projected displacements for a set of offline training cases. It is adapted for
%   **nonlinear finite element simulations in the large deformation regime** using the **CECM/NECM methodology**.
%
%   The stresses are evaluated on selected Gauss points (`INFO_RVE.GaussINDEX_stress`) and displacements on
%   selected DOFs (`INFO_RVE.DOFS_globNUM`). Both projected and exact stress snapshots are collected for later
%   ROM compression, error analysis, and model calibration.
%
% USAGE:
%   [SNAPreactFORCESproj, SNAPstressSTWOproj, SNAPstressPonePROJ, NAME_BASE, ...
%    SNAPstressPONEproj_NONLIN, SNAPstressPONEproj_LINEAR] = ...
%     GetStressesAndReactForces_bubNECMlarg(NAME_BASE, DATAcommon, NAMEsnap_base, DATAoffline, ...
%                                           INFO_RVE, BasisDEFdisp, DATA, OPERFE, ...
%                                           MATPRO, Fbody, Ftrac, OTHER_output, PhiRB)
%
% INPUTS:
%   - NAME_BASE         : Base name of the parametric study (used for labeling).
%   - DATAcommon        : Structure containing input parameters for all training cases.
%   - NAMEsnap_base     : String base for snapshot folder paths (per training sample).
%   - DATAoffline       : Structure with optional additional test cases and settings.
%   - INFO_RVE          : Structure specifying Gauss point indices and DOFs of interest for projection.
%   - BasisDEFdisp      : Structure with reduced displacement basis (e.g., `PhiDEFbs`, `PhiDEFcomp`).
%   - DATA              : Main data structure (simulation control, material flags, etc.).
%   - OPERFE            : Finite element operators and mappings.
%   - MATPRO            : Structure with material parameters.
%   - Fbody, Ftrac      : Cell arrays with body and traction force definitions for each case.
%   - OTHER_output      : Structure containing problem-specific output, e.g., precomputed FE variables.
%   - PhiRB             : Matrix with rigid body modes used in the displacement basis.
%
% OUTPUTS:
%   - SNAPreactFORCESproj        : Cell array of reaction force snapshots from projected displacements.
%   - SNAPstressSTWOproj         : Cell array of projected second Piola–Kirchhoff stress snapshots.
%   - SNAPstressPonePROJ         : (legacy/unused) Placeholder for PK1 stress projection.
%   - NAME_BASE                  : Passed-through input name for convenience.
%   - SNAPstressPONEproj_NONLIN  : Cell array of nonlinear (inelastic) PK1 projected stresses.
%   - SNAPstressPONEproj_LINEAR  : Cell array of linear (elastic) PK1 projected stresses.
%
% FUNCTIONALITY:
%   - Concatenates all displacement bases and builds projection operator via SVD.
%   - Loops over all training and additional test cases, loading snapshots and computing:
%       * Projected displacements via reduced basis
%       * Corresponding projected PK1 and PK2 stresses
%       * Reaction forces at selected DOFs
%   - Calls `MultiSnapStressFromDispNECMlarg` to evaluate nonlinear internal forces and stresses from displacement input.
%   - Computes and logs relative projection error of PK2 stresses for each case.
%
% FEATURES:
%   - Compatible with hybrid ROM strategies using both elastic and inelastic stress components.
%   - Supports additional offline test cases for validation or enrichment (via `DATAoffline.AdditionalTests`).
%   - Modular separation of projection (basis compression) and stress computation.
%   - Can be used for training of hyper-reduction, snapshot POD, or machine learning surrogate models.
%
% REFERENCES:
%   - Development based on:
%     /TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
%   - Theoretical context from:
%     /PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
%
% AUTHOR:
%   Joaquín A. Hernández, UPC – CIMNE
%   Date: 23-May-2024, Campus Nord, Barcelona
%   Comments by ChatGPT-4, 12-May-2025
%
% DEPENDENCIES:
%   - MultiSnapStressFromDispNECMlarg
%   - SVDT (for reduced basis projection)
%   - DefaultField
%   - INFO_SNAPSHOTS*.mat files for each case
%
% ---------------------------------------------------------------------------------------------------




% Copy of GetStressesAndReactForces_bubNECM.m
% Adaption to CECM just for nonlinear PK1 stresses, in large deformations
% JAHO, 23-May-2024, Barcelona, Campus Nord.  UPC
% SEe /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/108_EIFEM_metamat/04_EIFEM_NONL.mlx
if nargin == 0
    load('tmp1.mat')
end

% it is not clear that "PhiRB" should be introduced here, specially when it
% is  a large rotation problem
[BasisU,~,~] = SVDT([PhiRB,BasisDEFdisp.PhiDEFbs,BasisDEFdisp.PhiDEFcomp]) ;



% ------------------------------------------------------------------------------------------------------
% Now we have to compute the PK2 stresses produced by these strains, and check
% that all blocks are approximated to an accuracy theshold
% DATAoffline.\errorSTRESS
CASES = 1:size(DATAcommon.INPUTS_PARAMETERS,2) ;



VAR = [] ; VARint_n = [] ;
DATA = DefaultField(DATA,'ListFieldInternalVariables',[] ) ;





DATAoffline = DefaultField(DATAoffline,'AdditionalTests',[] );

if ~isempty(DATAoffline.AdditionalTests) ;
    CASES = 1:(size(DATAcommon.INPUTS_PARAMETERS,2) + length(DATAoffline.AdditionalTests)) ;
end

STRESS_PK2_error = zeros(length(CASES),1) ;

SNAPstressSTWOproj = cell(1,length(CASES)) ;
SNAPstressPONEproj_NONLIN = cell(1,length(CASES)) ; % "INELASTIC" STRESSES
SNAPstressPONEproj_LINEAR = cell(1,length(CASES)) ; % "ELASTIC" STRESSES

SNAPstressPonePROJ = cell(1,length(CASES)) ;
SNAPreactFORCESproj = cell(1,length(CASES)) ;

DATAoffline= DefaultField(DATAoffline,'SCALE_INFLUENCE_FACT_loc',1) ; 
SCALE_INFLUENCE_FACT_loc = DATAoffline.SCALE_INFLUENCE_FACT_loc ; 
for iproj = 1:length(CASES)
    % NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
    
    if iproj <= size(DATAcommon.INPUTS_PARAMETERS,2)
        NAME_FOLDER = [NAMEsnap_base,num2str(CASES(iproj))] ;
        ISADD = 0 ;
    else
        iloc = iproj-size(DATAcommon.INPUTS_PARAMETERS,2) ;
        NAME_FOLDER = DATAoffline.AdditionalTests{iloc} ;
        ISADD = 1;
    end
    
    NAME_INFO = [NAME_FOLDER,filesep,'INFO_SNAPSHOTS','.mat'] ;
    load(NAME_INFO,'INFO_SNAPSHOTS')
    
    DATA  = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA ; % 23-May-2024
    
    
    
    STORE_INFO = INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.STORE ;
    NAME_SNAP_loc = STORE_INFO.NAME_MATFILE_STORE ;
    SNAPstressSTWOproj_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SNAPstressPONEproj_LOC_NONL = cell(1,length(NAME_SNAP_loc)) ; % NONLINEAR STRESSES
    SNAPstressPONEproj_LOC_LINEAR = cell(1,length(NAME_SNAP_loc)) ; % LINEAR STRESSES
    
    SNAPstressSTWO_LOC  = cell(1,length(NAME_SNAP_loc)) ;
    
    
    
    SNAPstressPonePROJ_LOC = cell(1,length(NAME_SNAP_loc)) ;
    SELECTED_ENTRIES_GAUSS = INFO_RVE.GaussINDEX_stress ;
    SELECTED_ENTRIES_DOFS = INFO_RVE.DOFS_globNUM ;
    
    if ISADD ==0 && ~isempty(Fbody)
        Fbody_loc = Fbody{iproj} ;
        Ftrac_loc = Ftrac{iproj} ;
    else
        % Additional training test
        load(INFO_SNAPSHOTS.EMPLOYED_INPUTS.DATA.FE_VARIABLES_NAMEstore,'OTHER_output') ;
        Fbody_loc =OTHER_output.Fbody ;
        Ftrac_loc =OTHER_output.Ftrac ;
        if ~any(Fbody_loc.U)
            Fbody_loc.U = [] ;
            Fbody_loc.a = [] ;
        end
        
    end
    
    for iloc = 1:length(NAME_SNAP_loc)
        Nameloc = NAME_SNAP_loc{iloc} ;
        load(Nameloc,'SNAP_cluster') ;
        
        % 2ND PK STRESSES (Cauchy stresses for small strains )
        SNAPstressSTWO_LOC{iloc} = bsxfun(@times,SNAP_cluster.PK2STRESS.U(SELECTED_ENTRIES_GAUSS,:)',SNAP_cluster.PK2STRESS.S)' ;
        SNAPstressSTWO_LOC{iloc} = SNAPstressSTWO_LOC{iloc}*SNAP_cluster.PK2STRESS.V' ;
        
        % Projected displacements
        % -----------------------
        
        coeff = BasisU'*SNAP_cluster.DISP.U(SELECTED_ENTRIES_DOFS,:) ;
        coeff =bsxfun(@times,coeff',SNAP_cluster.DISP.S)' ;
        coeff = coeff*SNAP_cluster.DISP.V' ;
        % ------------------------------------------------
        d = BasisU*coeff; % Snapshot displacements (local)
        
        
        % Computation of stresses and reactive forces (residual) arising
        % from the projected stresses. Valid for both linear and  nonlinear
        % regimes
        [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR,SNAPreactFORCESproj{iproj},SNAPstressPONEproj_LOC_NONL,SNAPstressPONEproj_LOC_LINEAR] = ...
            MultiSnapStressFromDispNECMlarg...
            (VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
            SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj,INFO_RVE,Fbody_loc,Ftrac_loc,SNAPstressPONEproj_LOC_NONL,SNAPstressPONEproj_LOC_LINEAR) ;
        
    end
    SNAPstressSTWOproj_LOC = cell2mat(SNAPstressSTWOproj_LOC) ;   % Approximate
    SNAPstressPONEproj_LOC_NONL = cell2mat(SNAPstressPONEproj_LOC_NONL) ;   % Approximate
    SNAPstressPONEproj_LOC_LINEAR = cell2mat(SNAPstressPONEproj_LOC_LINEAR) ;   % Approximate
    
    SNAPstressSTWO_LOC = cell2mat(SNAPstressSTWO_LOC) ;   % exact (from FE training tests)
    
    % Recall that we may not want to study the whole range of snapshots ()
    %     nsnapshots = ceil(size(SNAPstressSTWO_LOC,2)*min(DATAoffline.proportionSNAPSHOTS(iproj),1)) ;
    %     SNAPstressSTWO_LOC = SNAPstressSTWO_LOC(:,1:nsnapshots) ;
    %     SNAPstressSTWOproj_LOC = SNAPstressSTWOproj_LOC(:,1:nsnapshots) ;
    %
    % Check if it meets the error criterion
    STRESS_PK2_error(iproj) = norm(SNAPstressSTWO_LOC-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;
    disp(['Project = ',num2str(iproj),'; ERROR stress PK2= ',num2str(STRESS_PK2_error(iproj))]);
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
    % PK2 stresses
    % ---------------
    %  [UU,SS,VV] =  RSVDT(SNAPreactFORCESproj{iproj}(:,1:nsnapshots)) ;
    %s  SNAPreactFORCESproj{iproj} = bsxfun(@times,UU',SS/SS(1))' ;
    
    
    SNAPstressSTWOproj{iproj} = SCALE_INFLUENCE_FACT_loc*SNAPstressSTWOproj_LOC ;
    SNAPstressPONEproj_NONLIN{iproj} = SCALE_INFLUENCE_FACT_loc*SNAPstressPONEproj_LOC_NONL ;
    SNAPstressPONEproj_LINEAR{iproj} = SCALE_INFLUENCE_FACT_loc*SNAPstressPONEproj_LOC_LINEAR ;
    
    
    
end



