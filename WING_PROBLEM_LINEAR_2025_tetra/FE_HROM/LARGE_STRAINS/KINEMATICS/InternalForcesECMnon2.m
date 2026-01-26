function  [Fint,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] = InternalForcesECMnon2(OPERFE,PoneST,PK2STRESS,DATA,VAR,tauNONallDER_q) ;
%--------------------------------------------------------------------------
% function [Fint, fINTredNON_mst_DOFl, wMSTeq, ResREDeqMST] = ...
%     InternalForcesECMnon2(OPERFE, PoneST, PK2STRESS, DATA, VAR, tauNONallDER_q)
%
% PURPOSE:
%   Computes the reduced internal force vector `Fint` in the context of a
%   nonlinear hyperreduced-order model (HROM) based on a manifold
%   representation. The contribution is separated into a "master" ECM-based
%   integration and a slave region approximation via a nonlinear map
%   (η_NON), both operating over the reduced coordinates of the system.
%
%   This function also computes auxiliary quantities used for tangent
%   system linearization:
%   - `wMSTeq`: equivalent integration weights (cf. Eq. (9.19))
%   - `ResREDeqMST`: residual vector defined with master weights (cf. Eq. (9.33))
%
% INPUTS:
% --------
%   OPERFE : struct
%       Contains offline precomputed quantities:
%       - Bst: stress projection operator at ECM points
%       - wSTs: integration weights for master points
%       - DOFl: indices of free degrees of freedom
%       - etaNON: nonlinear map for internal force extrapolation
%       - etaNONder: derivative of etaNON
%       - wRED_slv: integration weights for slave domain contribution
%
%   PoneST : [ngaus_master * ndim^2 × 1]
%       Vector of 1st Piola–Kirchhoff stress components at master ECM points.
%
%   PK2STRESS : unused if `PoneST` is non-empty; otherwise used for backup
%       Contains 2nd Piola–Kirchhoff stress components.
%
%   DATA : struct
%       Contains model and mesh definitions including spatial dimension `ndim`.
%
%   VAR : struct
%       Contains the current external force vector (`FEXT`) and other
%       quantities for tangent computations.
%
%   tauNONallDER_q : [n_modes_all × n_modes_manifold]
%       Derivative of manifold parametrization mapping, used to project
%       internal force densities onto the tangent space of the manifold.
%
% OUTPUTS:
% ---------
%   Fint : [n_modes_all × 1]
%       Reduced internal force vector including both master and extrapolated slave contributions.
%
%   fINTredNON_mst_DOFl : [n_modes_manifold × n_master_points]
%       Projected reduced internal force density at the master points,
%       restricted to free DOFs.
%
%   wMSTeq : [n_master_points × 1]
%       Equivalent integration weights that incorporate the extrapolated
%       contribution of the slave domain via η_NON'.
%
%   ResREDeqMST : [n_modes_all × 1]
%       Residual based on the equivalent integration at master points,
%       used in Jacobian/tangent computations (cf. Eq. (9.33)).
%
% REMARKS:
% --------
%   - This function implements Eq. (9.29) and (9.33) from the theory document.
%   - The nonlinearity of the reduced model is handled via a surrogate map η_NON.
%   - Slave contributions are indirectly accounted for using nonlinear
%     extrapolation based on the projected master internal force densities.
%   - If PoneST is empty, a fallback computation using PK2 stresses is provided
%     for compatibility/testing.
%
% DATE: 22-July-2025
% AUTHOR: Joaquín A. Hernández Ortega
%--------------------------------------------------------------------------

%
if nargin == 0
    load('tmp1.mat')
end


if ~isempty(PoneST)
    [Fint,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] = InternalForcesECMnon2LT(OPERFE,PoneST,DATA,VAR,tauNONallDER_q)  ;
else
    
    if  DATA.CECM_ONLY_FOR_NONLINEAR_STRESSES == 0
        error('Option not implemented yet, 29-August-2025')
        
        nF = DATA.MESH.nstrain ;
        for icomp = 1:nF
            icol = icomp:nF:length(PK2STRESS) ;
            PK2STRESS(icol,:) = PK2STRESS(icol,:).*OPERFE.wSTs;
        end
        Fint = OPERFE.Bst'*PK2STRESS ;
    else
        %         % Special implementation for EIFEM
        %         % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
        %         nF = DATA.MESH.nstrain ;
        %         for icomp = 1:nF
        %             icol = icomp:nF:length(PK2STRESS) ;
        %             PK2STRESS(icol,:) = VAR.PK2STRESS_incre(icol,:).*OPERFE.wSTs;
        %         end
        %         Fint = OPERFE.KstiffLINEAR*VAR.DISP + OPERFE.Bst'*PK2STRESS ;
        
        [Fint,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST,der_fINTredNON_slv_DOFl] = InternalForcesECMnon2SP(OPERFE,PK2STRESS,DATA,VAR,tauNONallDER_q)  ;
        
        
    end
end

