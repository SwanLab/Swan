function CECMoutput = CECM_mass_matrix(OPERFE,Phi,MATPRO,DATA,MESH,DATA_ECM)
%--------------------------------------------------------------------------
% function CECMoutput = CECM_mass_matrix(OPERFE, Phi, MATPRO, DATA, MESH, DATA_ECM)
%
% PURPOSE:
%   Constructs a reduced-order quadrature rule (via CECM) for efficiently 
%   integrating inertial (mass-related) forces in a nonlinear reduced-order model.
%
%   The goal is to select a small subset of integration points (empirical cubature)
%   that preserves the accuracy of mass matrix contributions, especially in
%   dynamics and model predictive simulations.
%
% INPUTS:
% -------
%   - OPERFE : structure with precomputed finite element operators, including:
%         · OPERFE.Nst     → shape functions at Gauss points (right-hand side)
%         · OPERFE.wSTs    → quadrature weights
%         · OPERFE.posgp_RHS → Gauss point coordinates (for mass integration)
%
%   - Phi : [ndof × nROM] reduced-order basis (displacement modes)
%
%   - MATPRO : material property structure, with field:
%         · MATPRO.dens → density vector (per element)
%
%   - DATA : global configuration with mesh and kinematics information
%
%   - MESH : mesh structure, including coordinates and Gauss point info
%
%   - DATA_ECM : structure configuring the empirical cubature algorithm
%
% OUTPUT:
% -------
%   - CECMoutput : structure containing:
%         · Selected integration points and weights (CECM)
%         · Projection error metrics
%         · Exact integral (full reference)
%         · DECM indices (if applicable)
%
% METHOD OVERVIEW:
% ----------------
% 1. Project the reduced basis `Phi` onto the finite element shape functions 
%    at Gauss points to form `NstRED`, representing reduced-order dynamics.
%
% 2. Construct the snapshot matrix for the integrand:
%        [ρ⋅NstRED, NstRED]
%    This allows the SVD to retain both inertial and geometric contributions.
%
% 3. Compute the reduced SVD basis `BasisForces` via truncated SVD (`SVDT`).
%
% 4. Use the routine `BasisA_InertialForces` to build the integrand matrix `A`,
%    which expresses the contribution of each reduced mode to the mass integral.
%
% 5. Apply the CECM algorithm (Continuous Empirical Cubature Method), which:
%      - Solves a sparse non-negative least squares problem
%      - Selects optimal Gauss points that reproduce the integrals of `A`
%
% 6. Optionally, evaluate the approximation error via `ErrorCalcLocal2023`.
%
% THEORY (CECM CONTEXT):
% ----------------------
%   This function implements the Continuous Empirical Cubature Method (CECM),
%   introduced in Hernández (2024), to generate accurate and sparse quadrature
%   rules for reduced-order modeling. Unlike classical EIM-based cubature
%   (which uses element-level interpolation), CECM works at the integrand level,
%   directly operating on pointwise values of internal force/inertial
%   contributions. This allows seamless integration with stress-based ROMs.
%
%   The integrand matrix `A` in this context satisfies:
%
%         ∫ ρ(x)⋅N(x)^T N(x) dx ≈ ∑_i w_i A(x_i)
%
%   where the nodes x_i and weights w_i are selected from the full integration set
%   using a sparsity-promoting algorithm.
%
% NOTES:
% ------
%   - This routine complements `CECM_internal_forces`, but focuses on mass terms.
%   - Designed for unstructured meshes and supports parallelization if needed.
%
% REFERENCES:
%   J. A. Hernández (2024), *A Continuous Empirical Cubature Method for
%   Nonlinear Model Order Reduction*, Computer Methods in Applied Mechanics and Engineering.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024
%--------------------------------------------------------------------------

if nargin ==0
    load('tmp1.mat')
end

% 
NstRED= OPERFE.Nst*Phi;  % M = int(Nst^T*W_RHS*Nst)
% Actually, NstRED is a matrix with as many rows as Gauss points time
% number of spatial dimensions, and with as many column and displacement
% modes. It represents the values of the nodal modes at RHS-Gauss points (right-hand side Gauss points)
% Density is defined at element level
dens = repmat(MATPRO.dens,1,DATA.MESH.ngaus_RHS*DATA.MESH.ndim)'  ;
dens = dens(:) ; 
NstREDdens = bsxfun(@times,NstRED,dens) ; 
% We have to include the density ---because our goal is to integrate exactly  inertial forces
ForcesAtGAussPoints = [NstREDdens,NstRED] ; 
DATASVD.ISRELATIVE = 1; 
[BasisForces,S,V] = SVDT(ForcesAtGAussPoints,1e-6,DATASVD) ; 
 


A = BasisA_InertialForces(NstRED,BasisForces,DATA)  ;
 

xNODES = MESH.COOR';
xFE = OPERFE.Nst*xNODES(:) ;

ndim = size(MESH.COOR,2) ;
xFE = reshape(xFE,ndim,[])' ;
%xFE = xFE(1:ndim:end,1) ;
%MESH.COOR = MESH.COOR(:,1)  ;
wFE = OPERFE.wSTs_RHS ;  % 30-May-2025




%delete('CECMoutputF.txt')
diary('CECMoutputF_b.txt')

% --------------------------------
% Evaluation integrands
% --------------------------------
 
DATAloc.ExactIntegral = A'*wFE ;
MESH.ngausE = DATA.MESH.ngaus_RHS ;
MESH.posgp = OPERFE.posgp_RHS ; % 10-May-2023
DATA_ECM.LabelPlot_integrand = 'BODYforces' ; 


% Continuous Empirical Cubature Method
[CECMoutput,DATA_AUX]= ContinuousECMgen2023(A,xFE,wFE,DATAloc,MESH,DATA_ECM) ;

CECMoutput.DECM_indexes_points   = DATA_AUX.indexPoints_DECM ; 

if DATA_ECM.UseDECMpoints ==0 

  ErrorCalcLocal2023(DATA,CECMoutput,DATA_AUX,A,DATAloc)  ;
end