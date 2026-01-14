function  [VAR,celastST,FgradST,detFgrad,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST] =  ResidualFromDisplacementsVARnonECMmd(OPERFE,VAR,MATPRO,DATA,VARint_n,tauNONallDER_q)
%--------------------------------------------------------------------------
% Adaptation of ResidualFromDisplacementsVARnonECM to cope with two-latent
% variable problems (elastoplasticy), and stresses decomposed into linear 
% and nonlinear 
% JAHO, 29th August 2025

% function [VAR, celastST, FgradST, detFgrad, fINTredNON_mst_DOFl, ...
%           wMSTeq, ResREDeqMST] = ...
%           ResidualFromDisplacementsVARnonECM(OPERFE, VAR, MATPRO, DATA, ...
%                                               VARint_n, tauNONallDER_q)
%
% PURPOSE:
%   Core routine for evaluating internal forces and residuals in nonlinear 
%   hyperreduced-order models (HROMs) with manifold-based representations.
%
%   The function:
%   - Evaluates stresses from current displacements using the constitutive law.
%   - Computes reduced internal forces using the master ECM integration points.
%   - Projects the internal forces onto the tangent space of the nonlinear manifold.
%   - Computes the residual vector (FINT - FEXT) in reduced coordinates.
%
%   It is a modified version of the classical residual evaluation routine, adapted
%   for use in nonlinear ROMs where the displacement field lies on a low-dimensional
%   nonlinear manifold parametrized by generalized coordinates.
%
% INPUTS:
% --------
%   OPERFE : struct
%       Contains precomputed finite element operators and ECM-related weights,
%       projection operators, nonlinear maps (`etaNON`), and DOF information.
%
%   VAR : struct
%       State variables at the current time step:
%       - DISP : nodal displacements
%       - FEXT : external force vector
%       - PK1STRESS, PK2STRESS : stress tensors to be updated
%
%   MATPRO : struct
%       Material parameters and constitutive law settings.
%
%   DATA : struct
%       Contains control flags for analysis, such as:
%       - INTERNAL_FORCES_USING_precomputed_Kstiff
%       - strain regime, mesh info, number of DOFs, etc.
%
%   VARint_n : struct
%       Internal variables at previous time step (needed for history-dependent
%       constitutive models, e.g., plasticity or viscoelasticity).
%
%   tauNONallDER_q : matrix [n_modes_all × n_modes_manifold]
%       Derivative of manifold parameterization (∂u/∂q) used to project
%       internal forces onto the tangent space of the reduced manifold.
%
% OUTPUTS:
% ---------
%   VAR : struct (updated)
%       - FINT : internal force vector in reduced space
%       - RESID : residual vector (FINT - FEXT)
%       - PK1STRESS, PK2STRESS : updated stress tensors
%
%   celastST : [nGauss × ndim^4] or Voigt-form elastic moduli tensor
%       Elastic stiffness tensor at all Gauss points, used for linearization.
%
%   FgradST : [nGauss × ndim^2]
%       Deformation gradient at each integration point (vectorized form).
%
%   detFgrad : [nGauss × 1]
%       Determinant of the deformation gradient.
%
%   fINTredNON_mst_DOFl : matrix
%       Reduced internal force density evaluated at master points, restricted
%       to free DOFs.
%
%   wMSTeq : vector
%       Equivalent weights for master integration points incorporating slave
%       region contributions via η_NON derivative.
%
%   ResREDeqMST : vector
%       Residual computed using only the master contribution with equivalent weights.
%       Used for tangent approximation.
%
% REMARKS:
% --------
%   - The function supports two modes:
%       (1) Full computation of stress/strain/deformation gradients, and
%       (2) A shortcut using a precomputed stiffness matrix (for linear models).
%
%   - For nonlinear HROMs, the function relies on `InternalForcesECMnon2`, which
%     implements the integration over master points and extrapolates the slave 
%     contribution via nonlinear maps (cf. Eq. 9.29–9.33 in theory).
%
%   - Essential for Newton–Raphson iterations in reduced-order nonlinear analysis.
%
% DATE: 22-July-2025
% AUTHOR: Joaquín A. Hernández Ortega
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp3.mat')
end



% Stresses from displacements
%[GLSTRAINS,PK2STRESS,celastST,FgradST,PoneST,detFgrad] = StressesFromDisplacements(OPERFE,d,MATPRO,DATA) ;
if DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==0
    [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsVAR(OPERFE,VAR,MATPRO,DATA,VARint_n) ;
else
    celastST = [] ; FgradST = [] ; detFgrad = [] ;
end


if isempty(VAR.PK2STRESS) && DATA.INTERNAL_FORCES_USING_precomputed_Kstiff ==0
    VAR.PoneST = [] ; Fint = [] ; VAR.RESID = [] ;
else
    OPERFE = DefaultField(OPERFE,'KinternalFORCES_given',[]) ;
    
  %  [VAR.FINT,fINTredNON_mst_DOFl,Bst_non_q] = InternalForcesECMnon(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA,VAR,tauNONallDER_q) ;
    % New linearization, 22-July-2025
        [VAR.FINT,fINTredNON_mst_DOFl,wMSTeq,ResREDeqMST] = InternalForcesECMnon2md(OPERFE,VAR.PK1STRESS,VAR.PK2STRESS,DATA,VAR,tauNONallDER_q,VAR.PK1STRESS_incre,...
            VAR.PK2STRESS_incre) ;

      VAR.FEXT = tauNONallDER_q'*VAR.FEXT ;
    %VAR.RESID = tauNONallDER_q'*VAR.RESID ;
    
    
    
    % 6.2. Residual
    VAR.RESID  = VAR.FINT- VAR.FEXT;
end


