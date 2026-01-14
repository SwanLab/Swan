function [SNAPreactFORCESproj,SNAPstressSTWOproj,SNAPstressPonePROJ,NAME_BASE,SNAPstressSTWOproj_INELAST,SNAPstressSTWOproj_ELAST]  =...
    GetStressesAndReactForces_bubNECM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
    INFO_RVE,BasisDEFdisp,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB)
 % =========================================================================
% FUNCTION: GetStressesAndReactForces_bubNECM
% =========================================================================
% PURPOSE:
% Computes projected stresses (2nd Piola-Kirchhoff and 1st PK), separated
% into elastic and inelastic parts, as well as projected reaction forces
% for each parameter configuration in the training set.
%
% This version is adapted to handle "bubble" modes (nonlinear complementary
% modes) in addition to rigid body and linear deformation modes.
%
% INPUTS:
% - NAME_BASE: Base name for tracking simulation files.
% - DATAcommon: Struct with shared problem data (parameter sets, etc.).
% - NAMEsnap_base: Base path to snapshot folders.
% - DATAoffline: Struct with settings for offline training and testing.
% - INFO_RVE: Information about the current RVE (e.g. Gauss points, DOFs).
% - BasisDEFdisp: Structure with basis fields for displacement approximation.
% - DATA: Problem data, including material and solver settings.
% - OPERFE: Finite element operator structure.
% - MATPRO: Material property data.
% - Fbody, Ftrac: Cell arrays of body and traction force structures per training case.
% - OTHER_output: Extra fields from full-order simulations.
% - PhiRB: Reduced basis (e.g. rigid body modes).
%
% OUTPUTS:
% - SNAPreactFORCESproj: Projected reaction forces per training sample.
% - SNAPstressSTWOproj: Projected 2nd Piola-Kirchhoff (PK2) stresses.
% - SNAPstressPonePROJ: Projected 1st Piola-Kirchhoff (PK1) stresses.
% - NAME_BASE: Updated base name (can reflect last read/write operation).
% - SNAPstressSTWOproj_INELAST: Inelastic part of PK2 stress projection.
% - SNAPstressSTWOproj_ELAST: Elastic part of PK2 stress projection.
%
% NOTES:
% - Uses a reduced displacement basis formed by combining PhiRB, linear
%   deformation modes (PhiDEFbs), and nonlinear complementary modes (PhiDEFcomp).
% - Calls `MultiSnapStressFromDispNECM` to compute projected stresses and
%   reaction forces from reduced displacements.
% - Computes a Frobenius-norm-based error metric to validate the accuracy of
%   stress projection.
% - Supports additional test cases via `DATAoffline.AdditionalTests`.
%
% RELATED:
% - Based on: GetStressesAndReactForces_1dom.m
% - Related testing script: /TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
% =========================================================================
% JAHO, updated 6-May-2025, comments by ChatGPT-4 
if nargin == 0
    load('tmp.mat')
end

 % BasisU is the basis matrix containing the rigid body modes, the linear
 % modes, and what we  call the complementary modes (nonlinear modes)
 [BasisU,~,~] = SVDT([PhiRB,BasisDEFdisp.PhiDEFbs,BasisDEFdisp.PhiDEFcomp]) ;
 
% ------------------------------------------------------------------------------------------------------
% Now we have to compute the PK2 stresses produced by these strains, and check
% that all blocks are approximated to an accuracy theshold
% DATAoffline.\errorSTRESS
CASES = 1:size(DATAcommon.INPUTS_PARAMETERS,2) ;

STRESS_PK2_error = zeros(length(CASES),1) ;

SNAPstressSTWOproj = cell(1,length(CASES)) ;
SNAPstressSTWOproj_INELAST = cell(1,length(CASES)) ; % "INELASTIC" STRESSES
SNAPstressSTWOproj_ELAST = cell(1,length(CASES)) ; % "ELASTIC" STRESSES

SNAPstressPonePROJ = cell(1,length(CASES)) ;
SNAPreactFORCESproj = cell(1,length(CASES)) ;

VAR = [] ; VARint_n = [] ;
DATA = DefaultField(DATA,'ListFieldInternalVariables',[] ) ;





DATAoffline = DefaultField(DATAoffline,'AdditionalTests',[] );

if ~isempty(DATAoffline.AdditionalTests) ;
    CASES = 1:(size(DATAcommon.INPUTS_PARAMETERS,2) + length(DATAoffline.AdditionalTests)) ;
end


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
    SNAPstressSTWOproj_LOC_INELAST = cell(1,length(NAME_SNAP_loc)) ; % INELASTIC STRESSES
    SNAPstressSTWOproj_LOC_ELAST = cell(1,length(NAME_SNAP_loc)) ; % ELASTIC STRESSES
    
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
        [SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,VAR,SNAPreactFORCESproj{iproj},SNAPstressSTWOproj_LOC_INELAST,SNAPstressSTWOproj_LOC_ELAST] = ...
            MultiSnapStressFromDispNECM...
            (VAR,d,DATA,VARint_n,OTHER_output,OPERFE,MATPRO,iloc,...
            SNAPstressPonePROJ_LOC,SNAPstressSTWOproj_LOC,iproj,INFO_RVE,Fbody_loc,Ftrac_loc,SNAPstressSTWOproj_LOC_INELAST,SNAPstressSTWOproj_LOC_ELAST) ;
        
    end
    SNAPstressSTWOproj_LOC = cell2mat(SNAPstressSTWOproj_LOC) ;   % Approximate
    SNAPstressSTWOproj_LOC_INELAST = cell2mat(SNAPstressSTWOproj_LOC_INELAST) ;   % Approximate
    SNAPstressSTWOproj_LOC_ELAST = cell2mat(SNAPstressSTWOproj_LOC_ELAST) ;   % Approximate
    
    SNAPstressSTWO_LOC = cell2mat(SNAPstressSTWO_LOC) ;   % exact (from FE training tests)
    
    
    % Check if it meets the error criterion
    STRESS_PK2_error(iproj) = norm(SNAPstressSTWO_LOC-SNAPstressSTWOproj_LOC,'fro')/norm(SNAPstressSTWO_LOC,'fro')  ;
    disp(['Project = ',num2str(iproj),'; ERROR stress PK2= ',num2str(STRESS_PK2_error(iproj))]);
    
    %%%% STORE SNAPSHOTS %%%%%%%% (COMPACT FORMAT, FOR DETERMINING AN ORTHOGONAL BASIS)
  
 
    if DATAoffline.BASIS_STRESSES_FROM_PROJECTED_DISPLACEMENTS ==1
        SNAPstressSTWOproj{iproj} = SNAPstressSTWOproj_LOC ;
        SNAPstressSTWOproj_INELAST{iproj} = SNAPstressSTWOproj_LOC_INELAST ;
        SNAPstressSTWOproj_ELAST{iproj} = SNAPstressSTWOproj_LOC_ELAST ;
    else
        error('Option not compatible')
        SNAPstressSTWOproj{iproj} = SNAPstressSTWO_LOC ;
    end
    
 
    
    
end
