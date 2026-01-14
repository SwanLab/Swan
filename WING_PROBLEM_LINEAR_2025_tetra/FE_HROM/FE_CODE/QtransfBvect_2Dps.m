function Be = QtransfBvect_2Dps(BeTILDE,ndim)
%--------------------------------------------------------------------------
% FUNCTION: QtransfBvect_2Dps
%
% PURPOSE:
%   Constructs the symmetric strain–displacement matrix `Be` for **plane strain**
%   conditions in 2D, including the out-of-plane strain component ε₃₃ = ε_zz.
%
%   This transformation maps nodal displacements to small strain components
%   in a 2D domain assuming generalized plane strain (i.e., with non-zero ε_z).
%
%   Given:
%       ε₁₁ = ∂u₁/∂x
%       ε₂₂ = ∂u₂/∂y
%       γ₁₂ = ∂u₁/∂y + ∂u₂/∂x
%       ε₃₃ = 0 (or prescribed separately)
%
%   The strain vector includes:
%       [ε₁₁, ε₂₂, γ₁₂, ε₃₃]ᵗ  → with ε₃₃ as a placeholder
%
%   Used primarily in J2 plasticity or multiphysics simulations that
%   require tracking ε_z in 2D sections (e.g., plane strain with 3D stress assumption).
%
% INPUTS:
%   - BeTILDE : Matrix of spatial derivatives of shape functions.
%               Size: (2 × nelem) × nnodeE
%   - ndim    : Number of spatial dimensions (must be 2)
%
% OUTPUT:
%   - Be      : Strain–displacement matrix of size (4 × nelem) × (2 × nnodeE)
%
% STRUCTURE:
%   Each block row of `Be` corresponds to a Gauss point and contains:
%     ε₁₁ → ∂N/∂x on u₁ DOFs
%     ε₂₂ → ∂N/∂y on u₂ DOFs
%     γ₁₂ → ∂N/∂y on u₁ + ∂N/∂x on u₂
%     ε₃₃ → (filled with zeros; treated as uncoupled or added later)
%
%   The implementation is a variant of `QtransfBvect`, with `nstrain = 4`
%   and ε_z treated explicitly.
%
% NOTES:
%   - The ε₃₃ component is often computed or imposed externally.
%   - The current implementation reserves a fourth row in each block for ε₃₃.
%   - `ndim` must be 2; if not, the behavior is undefined.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (jhortega@cimne.upc.edu)
%   CIMNE – Universitat Politècnica de Catalunya
%   Created: 8-Feb-2022
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp2.mat')
end
%Convertion from B-matrix for scalar fielts to B-matrix for vector fields
% Plane strain, last component = sigmaz, eZ 
%b rowcol = [1 2  6  3] ; 

nnodeE = size(BeTILDE,2); 
nelem = size(BeTILDE,1)/ndim ;

nstrain = 4 ; % Number of components  

Be = zeros(nstrain*nelem,nnodeE*ndim) ;
column1 = 1:2:(nnodeE*2-1) ;
column2 = 2:2:nnodeE*2 ;
ind1 = 1:nstrain:nstrain*nelem; ind2 = 2:nstrain:nstrain*nelem; ind3 = 3:nstrain:nstrain*nelem ;
jnd1 = 1:ndim:ndim*nelem; jnd2 = 2:ndim:ndim*nelem; 
Be(ind1,column1) = BeTILDE(jnd1,:) ;
Be(ind2,column2) = BeTILDE(jnd2,:) ;
Be(ind3,column1) = BeTILDE(jnd2,:) ;
Be(ind3,column2) = BeTILDE(jnd1,:) ;
