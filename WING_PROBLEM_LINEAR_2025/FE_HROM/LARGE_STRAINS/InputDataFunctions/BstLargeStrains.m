function [Bst_F,wSTs,ngaus,IDENTITY_F,posgp] = BstLargeStrains(MESH,nstrain)
%--------------------------------------------------------------------------
% FUNCTION: BstLargeStrains
%
% PURPOSE:
%   Constructs the global B-matrix `Bst_F` used for computing the
%   deformation gradient in the large strain regime, i.e.,
%
%       F(x) = I + Bst_F * d
%
%   where `d` is the global displacement vector, and `F` is the deformation
%   gradient evaluated at Gauss points. This B-matrix is tailored to
%   compute the deformation gradient from nodal displacements.
%
% INPUT:
%   - MESH: structure with finite element mesh information, including:
%       * COOR: nodal coordinates
%       * CN: element connectivity
%       * TypeElement: string indicating the type of element (e.g., 'Quadrilateral')
%       * posgp_given / weights_given: optional Gauss point positions/weights
%   - nstrain: number of strain components (e.g., 4 for 2D, 9 for 3D)
%
% OUTPUT:
%   - Bst_F: global matrix such that F = I + Bst_F*d
%   - wSTs : vector of Gauss weights times the Jacobian determinant
%   - ngaus: number of Gauss points per element
%   - IDENTITY_F: vectorized identity tensor (repeated at each Gauss point)
%   - posgp: positions of Gauss points (in parent domain)
%
% WHAT THIS FUNCTION DOES:
%   1. Constructs local matrices `Belem_F` used to compute the gradient of shape
%      functions at Gauss points.
%   2. Assembles them globally into `Bst_F`.
%   3. Builds a vectorized identity tensor `IDENTITY_F` for later use in
%      constructing F = I + B*d.
%
% REMARKS:
%   - This function supports both 2D and 3D meshes.
%   - The deformation gradient is computed at Gauss points, not at nodes.
%   - `Bst_F` maps global nodal displacements into gradients of displacements
%     (i.e., into the J = F - I tensor).
%
% REFERENCES:
%   See documentation file: DOCS/02_DeformationGradient_LagrangeStr.pdf
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC - CIMNE, Barcelona
%   Date: 22-Nov-2020
%--------------------------------------------------------------------------



% See DOCS/02_DeformationGradient_LagrangeStr.pdf
% % Bmatrix. Matrix such that DeformationGradient = Identity + Bst_F*d  
% wSTs is the vector of Gauss weights time the Jacobians (Jacobian of the transformation 
% from parent domain to physical domain )
% JAHO, 22-Nov-2020
if nargin ==0 
    load('tmp.mat')
end
% 
DATALOC.DEFORMATION_GRADIENT  = 1 ; 
MESH = DefaultField(MESH,'posgp_given',[]) ; 
MESH = DefaultField(MESH,'weights_given',[]) ; 

DATALOC.posgp_given = MESH.posgp_given ; 
DATALOC.weights_given = MESH.weights_given ; 

[ Belem_F, wSTs, wST, XeALL, posgp] = ComputeBelemALL(MESH.COOR,MESH.CN,MESH.TypeElement,nstrain,DATALOC) ;
% Computing the deformation gradient as a vector 
% F = J + I , where J = Bst_F*d
[nelem,nnodeE ]= size(MESH.CN) ; 
ngaus = length(wSTs)/nelem ;
disp('Assembly of Bst (stacked B-matrix)...')
[nnode,ndim ]= size(MESH.COOR) ;
if ndim == 2
    nstrainF = 4;
    IDENTITY = [1;1;0;0] ;
elseif ndim ==3 
    nstrainF = 9 ; 
    IDENTITY = [1;1;1;0;0;0;0;0;0] ; 
end
Bst_F = AssemblyBGlobal(Belem_F,nstrainF,nelem,nnodeE,ndim,ngaus,MESH.CN,nnode) ;


IDENTITY_F  = repmat(IDENTITY,1,nelem*ngaus) ; 
IDENTITY_F = IDENTITY_F(:) ; 


