function  [AsnapSTRESS_PK1,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel] ...
    = RetrievingSTRESSESprojected(DATAoffline,NAME_BASE,DATAcommon,NAMEsnap_base,INFO_RVE,BasisUdeform,DATA,...
    OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB)
%--------------------------------------------------------------------------
% FUNCTION: RetrievingSTRESSESprojected
%
% PURPOSE:
%   Retrieves and optionally decomposes the stress snapshots (e.g., PK1)
%   from training simulations, projected onto the reduced deformation basis.
%   The function handles both single-domain and multi-domain settings, and
%   distinguishes between linear and nonlinear contributions when required.
%
%   It supports two use cases:
%     1) Standard ROM training with one domain (default mode)
%     2) Multi-domain training (when auxiliary RVE data is provided)
%
%   Depending on `DATAoffline.CECM_ONLY_FOR_NONLINEAR_STRESSES`, the function
%   separates stress snapshots into elastic (`AsnapSTRESSel`) and inelastic
%   (`AsnapSTRESSinel`) parts for use in hyperreduction and basis truncation.
%
% INPUTS:
%   - DATAoffline    : Offline configuration (e.g., stress filtering flags)
%   - NAME_BASE      : Base string for naming and accessing snapshot files
%   - DATAcommon     : Configuration shared across training domains
%   - NAMEsnap_base  : Full path to snapshot storage folder
%   - INFO_RVE       : Metadata from the displacement training stage
%   - BasisUdeform   : Deformational modes computed from displacement data
%   - DATA           : General model and simulation data
%   - OPERFE         : Finite element operator structure
%   - MATPRO         : Material property data (e.g., constitutive tensors)
%   - Fbody/Ftrac    : External body forces and tractions
%   - OTHER_output   : Auxiliary outputs from displacement stage
%   - PhiRB          : Rigid body displacement modes
%
% OUTPUTS:
%   - AsnapSTRESS_PK1 : Projected PK1 stress snapshots (full or combined)
%   - NAME_BASE       : Possibly updated snapshot name base
%   - AsnapSTRESSinel : Inelastic (nonlinear) stress components (if applicable)
%   - AsnapSTRESSel   : Elastic stress components (if applicable)
%
% NOTES:
%   - If `DATA.INFO_RVEaux` is empty → single-domain case is assumed.
%   - If `CECM_ONLY_FOR_NONLINEAR_STRESSES == 0`, only full PK1 stresses
%     are returned (no decomposition).
%   - Stress decomposition methods differ for small strain vs. large strain setups.
%
% CALLED FUNCTIONS:
%   - GetStressesAndReactForces_bub
%   - GetStressesAndReactForces_bubNECM / bubNECMlarg
%   - RetrievingSTRESSESprojectedMULTI (for multi-domain setup)
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC-CIMNE)
%   First version: 7-Jan-2023 — Updated: 6-May-2025
%   Comments by ChatGPT4, 1-June-2025
%--------------------------------------------------------------------------


if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'INFO_RVEaux',[]) ; % Information "auxiliar" domains

if isempty(DATA.INFO_RVEaux)
    % Just one domain (option before version 6-May-2025)
    DATAoffline = DefaultField(DATAoffline,'CECM_ONLY_FOR_NONLINEAR_STRESSES',0) ; % 6-Jan-2023,  see
    % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
    if DATAoffline.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
        AsnapSTRESSinel = [] ;
        AsnapSTRESSel =[] ;
        [~ ,~ ,AsnapSTRESS_PK1,NAME_BASE]  = GetStressesAndReactForces_bub(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
            INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
    else
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        % Separate treatment nonlinear stresses
        if  DATA.SMALL_STRAIN_KINEMATICS ==1  && DATA.NO_USE_Deformation_gradient_in_Small_Strains == 1
            [~ ,~ ,AsnapSTRESS_PK1,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel]  =...
                GetStressesAndReactForces_bubNECM(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
                INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
        else
            [~ ,~ ,AsnapSTRESS_PK1,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel] ...
                = GetStressesAndReactForces_bubNECMlarg(NAME_BASE,DATAcommon,NAMEsnap_base,DATAoffline,...
                INFO_RVE,BasisUdeform,DATA,OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
        end
    end
    
else
    
    
    [AsnapSTRESS_PK1,NAME_BASE,AsnapSTRESSinel,AsnapSTRESSel] ...
        = RetrievingSTRESSESprojectedMULTI(DATAoffline,NAME_BASE,DATAcommon,NAMEsnap_base,INFO_RVE,BasisUdeform,DATA,...
        OPERFE,MATPRO,Fbody,Ftrac,OTHER_output,PhiRB) ;
    
    
end