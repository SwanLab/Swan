function  K =  KstifflDynamicPart(VAR,OPERFE,TIMEPARAM,DeltaT,K)
%--------------------------------------------------------------------------
% KstifflDynamicPart
% ------------------
% PURPOSE
%   Add **dynamic contributions** (effective mass and damping) to the tangent
%   matrix for Newmark–Bossak time integration:
%       K_eff  =  K  +  ((1-α)/(β Δt^2)) M  +  (γ/(β Δt)) C
%   where M ≡ OPERFE.M, C ≡ OPERFE.Ddamp (optional), and {α,β,γ} are
%   the Bossak/Newmark parameters stored in TIMEPARAM.
%
% INPUTS
%   VAR        : State container (not used here; passed for interface symmetry).
%   OPERFE     : FE/HROM operators:
%                • M      – (reduced) mass matrix
%                • Ddamp  – (optional) Rayleigh/viscous damping matrix; [] if none
%   TIMEPARAM  : Time integration parameters
%                • TYPE   – 1 ⇒ Newmark–Bossak
%                • a      – Bossak α (often ≤ 0)
%                • BETA   – Newmark β
%                • GAMMA  – Newmark γ
%   DeltaT     : Time step size (Δt)
%   K          : Current (elastic/tangent) stiffness matrix to augment
%
% OUTPUT
%   K          : Effective tangent including dynamic contributions
%
% FORMULA (TYPE == 1, Newmark–Bossak)
%   K ← K + ((1 - a)/(BETA * Δt^2)) * M  +  (GAMMA/(BETA * Δt)) * C,  with C=OPERFE.Ddamp if provided.
%
% NOTES & TIPS
%   • Units: ensure M and C are consistent with the discretization and Δt.
%   • Symmetry: if M and C are symmetric and K is symmetric, K_eff remains symmetric.
%   • Damping: if OPERFE.Ddamp is empty, only the mass term is added.
%   • Parameters: common unconditionally-stable choices are γ ≈ 0.5 - a and
%     β ≈ (1 - a)^2 / 4 (Bossak). Ensure TIMEPARAM matches your scheme.
%   • TYPE ≠ 1 is not implemented and will raise an error.
%
% SEE ALSO
%   Newmark/Bossak scheme setup, mass and damping assembly routines.
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp2.mat')
end



if TIMEPARAM.TYPE==1
    % Newmark-Bossak
%   \begin{equation}
%  \Ke{} =  \Ke{} + {\Par{\dfrac{1-\aBOSSAK}{\beta \Delta t_{n+1}^2} \M  +  \dfrac{\gamma}{\beta \Delta t_{n+1}} \Ddamp }  } 
% \end{equation}
    
    % Contribution displacement-independent damp./iner. external forces
     
    % Displacement-dependent term
    % ------------------------------
    K = K +  (1-TIMEPARAM.a)/(TIMEPARAM.BETA*DeltaT^2)*OPERFE.M ;
    if ~isempty(OPERFE.Ddamp)
         K =  K + TIMEPARAM.GAMMA/(TIMEPARAM.BETA*DeltaT)*OPERFE.Ddamp ;
    end
     
    
else
    error('Option not implemented')
    
    
end