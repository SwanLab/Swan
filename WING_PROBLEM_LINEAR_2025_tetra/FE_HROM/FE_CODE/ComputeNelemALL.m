function[ Nelem,posgp ] = ComputeNelemALL(TypeElement,nnodeE,ndim,nelem,TypeIntegrand) 
%--------------------------------------------------------------------------
% FUNCTION: ComputeNelemALL
%
% PURPOSE:
%   This function computes the matrix `Nelem`, which stacks the shape
%   function matrices (for vector-valued fields) for all elements at all
%   Gauss points. It is used primarily for assembling quantities such as
%   body forces, RHS vectors, or consistent mass matrices.
%
%   The returned matrix is independent of the geometry of the mesh because
%   it only depends on shape functions in the parent (reference) domain.
%
%   The matrix structure:
%       Nelem ∈ ℝ^{(nelem × ngaus × ndim) × (nnodeE × ndim)}
%
% INPUTS:
%   - TypeElement : (string) Type of finite element (e.g., 'Quadrilateral')
%   - nnodeE      : Number of nodes per element
%   - ndim        : Number of spatial dimensions (2 or 3)
%   - nelem       : Number of finite elements in the domain
%   - TypeIntegrand (optional): Type of integration ('RHS', 'mass', etc.).
%                               If not provided, defaults to 'RHS'.
%
% OUTPUTS:
%   - Nelem : Block matrix stacking the shape functions for all elements.
%             Each block is identical, because shape functions are the same
%             in the reference domain for all elements.
%   - posgp : Positions of Gauss points in the parent element.
%
% WHAT THIS FUNCTION DOES:
%   1. It first computes the shape functions and Gauss points in the
%      reference domain using `ComputeElementShapeFun`.
%   2. It builds the local shape function matrix `Ne` for a vector-valued
%      field at each Gauss point, using `StransfN`.
%   3. Since these shape functions are invariant across elements in the
%      reference configuration, it constructs `Nelem` by replicating `Ne`
%      for each element using `repmat`.
%
% REMARKS:
%   - This function does **not** depend on the actual geometry (COOR or CN).
%   - The output matrix is typically used in routines for assembling
%     consistent mass matrices or evaluating distributed forces.
%   - This is a **vectorized** version for efficiency.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   CIMNE – Universitat Politècnica de Catalunya
%   Date: 26-Oct-2015
%   Comments by ChatGPT4, 12-May-2025
%--------------------------------------------------------------------------

if nargin == 4
    TypeIntegrand = 'RHS';
end
%%%%
% This subroutine   returns
%  1) the matrix Nelem (ndim*nelem*ngaus x ndim*nnodeE)  consisting of all element shape function matrices
% Inputs:   COOR: Coordinate matrix (nnode x ndim), % CN: Connectivity matrix (nelem x nnodeE),
% TypeElement: Type of finite element (quadrilateral,...),
%% Vectorized version
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 26-Oct-2015
%dbstop('10')
if nargin == 0
    load('tmp2.mat')
end
%nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;
% Shape function routines (for calculating shape functions and derivatives)


[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ;
ngaus = length(weig) ;
% Shape functions for vector  fields 
Ne = zeros(ngaus*ndim,nnodeE*ndim) ; 
for g = 1:ngaus
    indI = (g-1)*ndim+1; indF = g*ndim ; 
    Ne(indI:indF,:) = StransfN(shapef(g,:),ndim) ;
end
% Shape functions are INDEPENDENT of the physical coordinates
%                     ---------------------------------------
% Therefore, the sought-after Nelem matrix can be computed by making  nelem tiling
% copies of Ne
Nelem = repmat(Ne,nelem,1); 



% Initialization
% Nelem = zeros(ndim*nelem*ngaus,nnodeE*ndim) ;
% % Let us define a matrix ROWSgauss such that   Belem(ROWSgauss(g,:),:) returns 
% % the N-matrices of the g-th points of all elements
% indREF = 1:nstrain*ngaus*nelem ;
% ROWSgauss = reshape(indREF,nstrain,nelem*ngaus) ; 
% weigREP = repmat(weig',nelem,1)  ; % nelem x 1 tiling copies of weig  
% for  g = 1:ngaus
%     % Matrix of derivatives for Gauss point "g"
%     BeXi = dershapef(:,:,g) ;
%     % Jacobian Matrix for the g-th G. point of all elements %
%     JeALL = XeALL*BeXi' ;
%     %%%%%%%%%
%     % JAcobian
%     detJeALL= determinantVECTORIZE(JeALL,ndim) ;
%     % Matrix of derivatives with respect to physical coordinates
%     inv_JeTall = inverseTRANSvectorize(JeALL,ndim,detJeALL) ;
%     BeTILDEall = inv_JeTall*BeXi ;
%     % Matrix of symmetric gradient
%     % dbstop('18')
%     BeALL = QtransfBvect(BeTILDEall,ndim) ; % Transformation from B-matrix for scalar fields to B-matrix for vector fields
%     % Assigning BeALL ( g-th Gauss point) to the global matrix Belem    
%     ROWSglo =  ROWSgauss(:,g:ngaus:ngaus*nelem);   
%     ROWSglo = ROWSglo(:) ;      
%     Belem(ROWSglo,:) = BeALL ;   
%     % Weight vectors
%     % --------------
%     [wST,wSTs] = DetermineWeightsST(detJeALL,weigREP,ngaus,g,nelem,wST,wSTs,nstrain,ROWSglo);
%     
% end
