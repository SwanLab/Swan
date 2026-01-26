function  [AsnapSTRESS_PK1_all,NAME_BASE,AsnapSTRESSinel_all,AsnapSTRESSel_all] ...
    = RetrievingSTRESSESprojectedMULTI(DATAoffline,NAME_BASE,DATAcommon,NAMEsnap_base,INFO_RVEref,BasisUdeform,DATA,...
    OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB)
% =========================================================================
% FUNCTION: RetrievingSTRESSESprojectedMULTI
% =========================================================================
% PURPOSE:
% This function computes and retrieves the projected stress fields (e.g.,
% 1st Piola-Kirchhoff stresses, elastic and inelastic stress components)
% for multiple finite element (FE) subdomains (reference and auxiliary).
%
% It supports both standard and nonlinear stress recovery depending on the configuration flags in
% DATAoffline and DATA. For nonlinear problems, it distinguishes between
% small strain and large rotation regimes.
%
% INPUTS:
% - DATAoffline: Structure containing offline settings for reduced order modeling.
% - NAME_BASE: Base name for saving/loading snapshots and data.
% - DATAcommon: General shared data for the FE problem.
% - NAMEsnap_base: Base filename for the snapshot files.
% - INFO_RVEref: Struct containing reference RVE domain information (nodes, DOFs, elements, Gauss points).
% - BasisUdeform: Deformation basis for the projection.
% - DATA: General data struct, including conversion indices for auxiliary domains.
% - OPERFE: Operator data for FE computations.
% - MATPRO: Material properties.
% - Fbody, Ftrac: External body and traction forces.
% - OTHER_output: Additional output controls and configuration.
% - PhiRB: Reduced basis for ROM projection.
%
% OUTPUTS:
% - AsnapSTRESS_PK1: Cell array of projected 1st Piola-Kirchhoff stresses for each subdomain.
% - NAME_BASE: Updated base name, if needed.
% - AsnapSTRESSinel: Cell array of inelastic stress components (if ECM for nonlinearity is enabled).
% - AsnapSTRESSel: Cell array of elastic stress components (if ECM for nonlinearity is enabled).
%
% NOTES:
% - For auxiliary domains, appropriate index mapping to the reference domain is
%   performed using DATA.ConversionINDEX_AUXDOM.
% - Functionality depends on flags such as CECM_ONLY_FOR_NONLINEAR_STRESSES
%   and SMALL_STRAIN_KINEMATICS to choose between standard or nonlinear
%   treatment (GetStressesAndReactForces_bub or GetStressesAndReactForces_bubNECM).
% - See associated test script:
%   /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
%
% =========================================================================
% JAHO, 6-May-2025, Terrassa UPC, comments by CHATGPT 4
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
if nargin == 0
    load('tmp1.mat')
end

% These are the indices that allows one to extract information of the FE
% results of the entire training domain

% ----------------------------------
% REFERENCE DOMAIN
% ----------------------------------
% INFO_RVEref
%
% INFO_RVEref =
%
%   struct with fields:
%
%            NODES_globNUM: [1548×1 double]
%             DOFS_globNUM: [3096×1 double]
%         BODYELEM_globNUM: [344×1 double]
%     GaussINDEX_scalarSTR: [3096×1 double]
%       GaussINDEX_stress: [12384×1 double]

% FOR THE "AUXILIAR" DOMAINS
%INFO_RVEaux = DATA.INFO_RVEaux;
%  INFO_RVEaux{1}
%
% ans =
%
%   struct with fields:
%
%            NODES_globNUM: [1548×1 double]
%             DOFS_globNUM: [3096×1 double]
%         BODYELEM_globNUM: [344×1 double]
%     GaussINDEX_scalarSTR: [3096×1 double]
%        GaussINDEX_stress: [12384×1 double]

% HOEWEVER, THE AUXILIAR DOMAINS NEED TO BE PUT INTO CORRESPONDENCE TO THE
% REFERENCE DOMAIN. THE CONVERSION DATA IS STORED IN
%
% DATA.ConversionINDEX_AUXDOM
%
% ans =
%
%   struct with fields:
%
%                NodalDOFS: {[3096×1 double]}
%                 ELEMENTS: {[344×1 double]}
%     GaussINDEX_scalarSTR: {[3096×1 double]}
%        GaussINDEX_stress: {[12384×1 double]}

nAUXdom  =length(DATA.INFO_RVEaux)  ;
nDOMtotal = nAUXdom ;
AsnapSTRESS_PK1 = cell(1,nDOMtotal) ;
AsnapSTRESSinel = cell(1,nDOMtotal) ;
AsnapSTRESSel = cell(1,nDOMtotal) ;
% This is for scaling the influence of the resulting stress snapshots 
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/ScalingSnapshotsDISP.m
SCALE_INFLUENCE_FACT = DATAoffline.INFLUENCE_SNAPSHOTS_DISP.FIRST_DOMAIN;
SCALE_INFLUENCE_FACT = [SCALE_INFLUENCE_FACT;DATAoffline.INFLUENCE_SNAPSHOTS_DISP.AUXDOMAINS(:)];

 
for    idom = 1:(nAUXdom+1)
    
    DATAoffline.SCALE_INFLUENCE_FACT_loc = SCALE_INFLUENCE_FACT(idom) ; 
    
    if idom ==1
        % REFERENCE DOMAIN
        fprintf('----------------------------------------------\n');
        fprintf('Stress error for the reference domain\n');
        fprintf('----------------------------------------------\n');
        INFO_RVE = INFO_RVEref ;
    else
        
        iaux = idom-1 ;
        fprintf('----------------------------------------------\n');
        fprintf('Stress error for the auxiliar domain = %d\n', iaux);
        fprintf('----------------------------------------------\n');
        
        INFO_RVE = DATA.INFO_RVEaux{iaux} ;
        % All we need to change is
        %  SELECTED_ENTRIES_GAUSS = INFO_RVE.GaussINDEX_stress ;
        % SELECTED_ENTRIES_DOFS = INFO_RVE.DOFS_globNUM ;
        % According to what it was done in
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/DOMAINdecompHROM/GENERAL_BUBBLE/GetDisplacAndMesh_sevDOM.m
        % for nodal DOFs
        %   SELECTED_DOFS = INFO_RVEaux{iaux}.DOFS_globNUM ;
        %    SELECTED_DOFS = SELECTED_DOFS(ConversionINDEX_AUXDOM.NodalDOFS{iaux}) ;
        %    for iloc = 1:length(NAME_SNAP_loc)
        %                 Nameloc = NAME_SNAP_loc{iloc} ;
        %                 load(Nameloc,'SNAP_cluster') ;
        %                 % Rather than reconstructing as U*S*V', we only make U*S (weighted left singular vectors)
        %                 DISP_LOC{iloc} = bsxfun(@times,SNAP_cluster.DISP.U(SELECTED_DOFS,:)',SNAP_cluster.DISP.S)' ;
        
        INFO_RVE.DOFS_globNUM  = INFO_RVE.DOFS_globNUM(DATA.ConversionINDEX_AUXDOM.NodalDOFS{iaux}) ;
        INFO_RVE.GaussINDEX_stress  = INFO_RVE.GaussINDEX_stress(DATA.ConversionINDEX_AUXDOM.GaussINDEX_stress{iaux}) ;
        
        
    end
    
    
    if DATAoffline.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
        AsnapSTRESSinel = [] ;
        AsnapSTRESSel =[] ;
        warning('This function has not been updated properly...Follow the steps given in GetStressesAndReactForces_bubNECMlarg')
        [~ ,~ ,AsnapSTRESS_PK1{idom},NAME_BASE]  = GetStressesAndReactForces_bub(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
            INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
    else
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        % Separate treatment nonlinear stresses
        if  DATA.SMALL_STRAIN_KINEMATICS ==1  && DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1
           warning('This function has not been updated properly (SCALE_INFLUENCE_FACT)...Follow the steps given in GetStressesAndReactForces_bubNECMlarg')   
            [~ ,~ ,AsnapSTRESS_PK1{idom},NAME_BASE,AsnapSTRESSinel{idom},AsnapSTRESSel{idom}]  =...
                GetStressesAndReactForces_bubNECM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
                INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
        else
            %
            [~ ,~ ,AsnapSTRESS_PK1{idom},NAME_BASE,AsnapSTRESSinel{idom},AsnapSTRESSel{idom}] ...
                = GetStressesAndReactForces_bubNECMlarg(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
                INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
        end
    end
    
    
end
disp('Arranging the stress snapshots')

numberPROJECTS= length(AsnapSTRESSel{1}) ;

AsnapSTRESS_PK1_all = cell(1,numberPROJECTS) ; 
AsnapSTRESSinel_all = cell(1,numberPROJECTS) ; 
AsnapSTRESSel_all  =  cell(1,numberPROJECTS) ; 


for idom = 1:length(AsnapSTRESS_PK1)
    for iproj = 1:numberPROJECTS
        AsnapSTRESS_PK1_all{iproj} = [AsnapSTRESS_PK1_all{iproj},AsnapSTRESS_PK1{idom}{iproj}] ; 
        AsnapSTRESSel_all{iproj} = [AsnapSTRESSel_all{iproj},AsnapSTRESSel{idom}{iproj}] ; 
         AsnapSTRESSinel_all{iproj} = [AsnapSTRESSinel_all{iproj},AsnapSTRESSinel{idom}{iproj}] ;
    end
end
 
 