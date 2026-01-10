function Be = QtransfBvect(BeTILDE,ndim)
%--------------------------------------------------------------------------
% FUNCTION: QtransfBvect
%
% PURPOSE:
%   Transforms the strain–displacement matrix `BeTILDE`, derived from scalar
%   derivatives of shape functions (∂N/∂x), into the full symmetric gradient
%   matrix `Be` used in linearized kinematics.
%
%   The resulting matrix `Be` maps nodal displacements to small strains.
%
%   This function implements the standard conversion:
%
%     ε = B * d
%
%   where:
%     - ε is the strain vector [ε₁₁, ε₂₂, γ₁₂]' in 2D or [ε₁₁, ε₂₂, ε₃₃, γ₂₃, γ₁₃, γ₁₂]' in 3D,
%     - d is the global nodal displacement vector,
%     - B is the symmetric gradient operator assembled here.
%
% INPUTS:
%   - BeTILDE : Matrix of shape function derivatives in physical coordinates.
%               Size: (ndim × nelem) × nnodeE
%   - ndim    : Spatial dimension (2 or 3)
%
% OUTPUT:
%   - Be      : Final symmetric B-matrix (nstrain × nelem) × (ndim × nnodeE)
%               Used in forming stiffness matrix K = Bᵗ * D * B
%
% DIMENSIONAL DETAILS:
%   - 2D case (ndim = 2): strain vector = [ε₁₁, ε₂₂, γ₁₂]
%     Assignments:
%       Be(1,:) = dN/dx
%       Be(2,:) = dN/dy
%       Be(3,:) = dN/dy on u₁ dofs + dN/dx on u₂ dofs (engineering γ₁₂)
%
%   - 3D case (ndim = 3): strain vector = [ε₁₁, ε₂₂, ε₃₃, γ₂₃, γ₁₃, γ₁₂]
%     Each row of Be corresponds to:
%       ε₁₁ → dN/dx on u₁
%       ε₂₂ → dN/dy on u₂
%       ε₃₃ → dN/dz on u₃
%       γ₂₃ → dN/dz on u₂ + dN/dy on u₃
%       γ₁₃ → dN/dz on u₁ + dN/dx on u₃
%       γ₁₂ → dN/dy on u₁ + dN/dx on u₂
%
%   - The resulting matrix is sparse and block-structured, with each strain
%     component assigned to the correct degrees of freedom.
%
% USAGE CONTEXT:
%   - Used in finite element assembly routines for small strain elasticity.
%   - Commonly called within `ComputeBelemALL` or equivalent routines.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (jhortega@cimne.upc.edu)
%   CIMNE – Universitat Politècnica de Catalunya
%   First version: 26-Oct-2015
%   Last revised: 22-Nov-2020
%   Comments by ChatGPT4, 15-May-2025
%--------------------------------------------------------------------------

%
if nargin == 0
    load('tmp2.mat')
end
%Convertion from B-matrix for scalar fielts to B-matrix for vector fields 

nnodeE = size(BeTILDE,2);
nelem = size(BeTILDE,1)/ndim ;
if ndim==2
    nstrain = 3 ;    Be = zeros(nstrain*nelem,nnodeE*ndim) ;
    column1 = 1:2:(nnodeE*2-1) ;
    column2 = 2:2:nnodeE*2 ;
    ind1 = 1:nstrain:nstrain*nelem; ind2 = 2:nstrain:nstrain*nelem; ind3 = 3:nstrain:nstrain*nelem ;
    jnd1 = 1:ndim:ndim*nelem; jnd2 = 2:ndim:ndim*nelem; jnd3 = 3:ndim:ndim*nelem ;
    Be(ind1,column1) = BeTILDE(jnd1,:) ;
    Be(ind2,column2) = BeTILDE(jnd2,:) ;
    Be(ind3,column1) = BeTILDE(jnd2,:) ;
    Be(ind3,column2) = BeTILDE(jnd1,:) ;
elseif ndim ==3
    nstrain = 6 ;    Be = zeros(nstrain*nelem,nnodeE*ndim) ;
    column1 = 1:3:(nnodeE*3-2) ;
    column2 = 2:3:(nnodeE*3-1) ;
    column3 = 3:3:nnodeE*3 ;
    ind1 = 1:nstrain:nstrain*nelem; ind2 = 2:nstrain:nstrain*nelem; ind3 = 3:nstrain:nstrain*nelem ;
    ind4 = 4:nstrain:nstrain*nelem; ind5 = 5:nstrain:nstrain*nelem; ind6 = 6:nstrain:nstrain*nelem ;
    jnd1 = 1:ndim:ndim*nelem; jnd2 = 2:ndim:ndim*nelem; jnd3 = 3:ndim:ndim*nelem ;
    jnd4 = 4:ndim:ndim*nelem; jnd5 = 5:ndim:ndim*nelem; jnd6 = 6:ndim:ndim*nelem ;
    Be(ind1,column1) = BeTILDE(jnd1,:) ;
    Be(ind2,column2) = BeTILDE(jnd2,:) ;
    Be(ind3,column3) = BeTILDE(jnd3,:) ;
    Be(ind4,column2) = BeTILDE(jnd3,:) ;
    Be(ind4,column3) = BeTILDE(jnd2,:) ;
    Be(ind5,column1) = BeTILDE(jnd3,:) ;
    Be(ind5,column3) = BeTILDE(jnd1,:) ;
    Be(ind6,column1) = BeTILDE(jnd2,:) ;
    Be(ind6,column2) = BeTILDE(jnd1,:) ;
else
    error('Incorrect option')
end