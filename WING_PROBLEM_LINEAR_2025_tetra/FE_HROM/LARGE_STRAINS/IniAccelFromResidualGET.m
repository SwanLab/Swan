function VAR = IniAccelFromResidualGET(DATA,DISP_CONDITIONS,OPERFE,Fbody,Ftrac,VAR,MATPRO)  
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Initialize accelerations at t = t₀ from the dynamic residual, enforcing
%   compatibility with prescribed (constrained) DOFs. The routine:
%     1) Infers constrained-DOF velocity at t₁ and acceleration at t₀ from
%        the prescribed displacements (finite difference on STEPS(1:2)).
%     2) Evaluates the static residual at t₀ from current displacements d,
%        internal forces (via constitutive update) and external loads.
%     3) Solves for free-DOF acceleration a₀ₗ from the block of the mass
%        matrix on free DOFs, including cross inertial and damping terms.
%     4) Assembles the full a₀ and updates state/diagnostic fields.
%
% INPUTS:
%   DATA                : problem data (uses DATA.STEPS for t₀, t₁).
%   DISP_CONDITIONS     : boundary-condition sets and data:
%                           • DOFl → free dof indices, DOFr → constrained dof indices
%                           • dR.U, dR.a(:,2) → prescribed displacement at step 2
%   OPERFE              : FE operators (uses M and optional Ddamp).
%   Fbody, Ftrac        : external load factors with F = U * a(:,istep).
%   VAR                 : current state (uses VAR.DISP = d, VAR.VEL = v).
%   MATPRO              : material data for internal force evaluation.
%
% OUTPUTS:
%   VAR                 : updated with
%                           • VAR.ACEL      → initial acceleration a₀
%                           • VAR.PK2STRESS → 2nd Piola–Kirchhoff stress at t₀
%                           • VAR.GLSTRAINS → Green–Lagrange strains at t₀
%                         plus energies via KineticAndStrainEnergyGet.
%
% METHOD (key steps and equations):
%   - Constrained kinematics (on DOFr):
%       d₁ʳ = dR.U * dR.a(:,2)
%       d₀ʳ = d(DOFr)
%       v₁ʳ = (d₁ʳ − d₀ʳ) / (t₁ − t₀)
%       a₀ʳ = (v₁ʳ − v(DOFr)) / (t₁ − t₀)
%
%   - External load at t₀ (istep=1):
%       FEXT = Fbody.U*a_body(:,1) + Ftrac.U*a_trac(:,1)
%
%   - Static residual on free DOFs at t₀:
%       RESID_static = (FINT(d₀) − FEXT)(DOFl)
%       (obtained via ResidualFromDisplacements, together with PK2STRESS/GLSTRAINS)
%
%   - Solve for free accelerations a₀ₗ:
%       Mll * a₀ₗ = − ( RESID_static + Mlr * a₀ʳ + Dll * v )
%       where Mll = M(DOFl,DOFl), Mlr = M(DOFl,DOFr),
%             Dll = Ddamp(DOFl,:) (if provided; otherwise 0)
%
%   - Assemble full acceleration:
%       a₀(DOFl) = a₀ₗ,   a₀(DOFr) = a₀ʳ
%
% ASSUMPTIONS / SANITY:
%   - DATA.STEPS has at least two entries; t₁ > t₀.
%   - DOFl and DOFr are disjoint and cover the relevant DOFs consistently.
%   - OPERFE.M is positive definite on DOFl; Ddamp may be empty.
%   - External loads Fbody/Ftrac provide coefficients for istep=1.
%
% NUMERICAL NOTES:
%   - Finite-difference estimates for v₁ʳ and a₀ʳ rely on Δt = t₁−t₀; very
%     small Δt amplifies round-off. Ensure STEPS spacing is reasonable.
%   - The constitutive update and FINT evaluation are encapsulated in
%     ResidualFromDisplacements.
%
% AUTHOR / HISTORY:
%   Comments clarification: 7-Nov-2025
%   JAHO — Joaquín A. Hernández — jhortega@cimne.upc.edu
% -------------------------------------------------------------------------

if nargin == 0
    load('tmp2.mat')
end

% Acceleration of the constrained DOFs
d = VAR.DISP ; 
v = VAR.VEL ; 

% Displacements 
d1r =  DISP_CONDITIONS.dR.U*DISP_CONDITIONS.dR.a(:,2) ; 
d0r = d(DISP_CONDITIONS.DOFr) ; 
 
t1 = DATA.STEPS(2) ; t0 = DATA.STEPS(1) ; 

% Velocity at time t=t_1 
v1r = (d1r-d0r)/(t1-t0) ; 

% Acceleration at time t_0 
a0r = (v1r-v(DISP_CONDITIONS.DOFr))/(t1-t0) ; 


%  RESIDUAL at t= 0
%  \begin{equation}
%  \label{eq:sdfas--}
%   \Res{}{\DOFl}(t_0) =  \M_{\DOFl \DOFl } \a^0_{\DOFl} + \M_{\DOFl \DOFr} \a^0_{\DOFr} +

% \Ddamp_{\DOFl} \vINI   +   \Fint{}{\DOFl}(\dINI)   -    \Fext{}{\DOFl}(t_{0})  = \zero
%  \end{equation}

% Mll a_0_l = -RESID_static -\M_{\DOFl \DOFr} \a^0_{\DOFr} -\Ddamp_{\DOFl} \vINI 

% Difference between external and internal forces 
% RESID_static = 
istep = 1; 
  FEXT = Fbody.U*Fbody.a(:,istep) + Ftrac.U*Ftrac.a(:,istep) ;
[GLSTRAINS,PK2STRESS,~,RESID_static,~,~,~,~] = ResidualFromDisplacements(OPERFE,d,MATPRO,DATA,FEXT)  ; 
 
  RESID_static= RESID_static(DISP_CONDITIONS.DOFl) ; 
  
% INERTIAL TERM % \M_{\DOFl \DOFr} \a^0_{\DOFr}
RESID_iner = OPERFE.M(DISP_CONDITIONS.DOFl,DISP_CONDITIONS.DOFr)*a0r ;  
% DAMPING TERM
if ~isempty(OPERFE.Ddamp)
    RESID_damp= OPERFE.Ddamp(DISP_CONDITIONS.DOFl,:)*v ;
else
    RESID_damp = 0 ; 
end

a0l = -OPERFE.M(DISP_CONDITIONS.DOFl,DISP_CONDITIONS.DOFl)\(RESID_static+RESID_iner+RESID_damp) ; 

a0 = zeros(size(d)) ; 
 a0(DISP_CONDITIONS.DOFl)  = a0l ; 
 a0(DISP_CONDITIONS.DOFr)  = a0r ; 

VAR.ACEL = a0 ; 
VAR.PK2STRESS = PK2STRESS ; 
VAR.GLSTRAINS = GLSTRAINS ; 


VAR = KineticAndStrainEnergyGet(VAR,OPERFE,DATA) ; 

