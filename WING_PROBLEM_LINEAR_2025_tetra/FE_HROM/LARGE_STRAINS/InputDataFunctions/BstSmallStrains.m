function [Bst_F,wSTs,ngaus,IDENTITY_F,posgp] = BstSmallStrains(MESH,nstrain)
%--------------------------------------------------------------------------
% FUNCTION: BstSmallStrains
%
% PURPOSE:
%   Constructs the global B-matrix `Bst_F` used to compute the linearized
%   small strain tensor from nodal displacements:
%
%       ε(x) = Bst_F * d
%
%   where `d` is the global displacement vector, and `ε(x)` is the small
%   strain tensor evaluated at Gauss points.
%
%   This function is the small strain counterpart to `BstLargeStrains`, and
%   is used when linearized kinematics are assumed.
%
% INPUT:
%   - MESH: structure with the finite element mesh information, including:
%       * COOR: nodal coordinates
%       * CN: element connectivity
%       * TypeElement: string with the type of finite element
%       * DATA.StrainStressWith4Components (optional): used to adapt nstrain
%   - nstrain: number of strain components (e.g., 3 for 2D, 6 for 3D)
%
% OUTPUT:
%   - Bst_F: global strain-displacement matrix
%   - wSTs : vector of Gauss weights times Jacobians
%   - ngaus: number of Gauss points per element
%   - IDENTITY_F: empty (not used in small strain formulation)
%   - posgp: positions of Gauss points in the reference (parent) domain
%
% WHAT THIS FUNCTION DOES:
%   1. Sets DEFORMATION_GRADIENT flag to 0 to disable nonlinear kinematics.
%   2. Calls `ComputeBelemALL` to compute local strain-displacement matrices.
%   3. Assembles local matrices into the global B-matrix `Bst_F`.
%   4. Returns weights, number of Gauss points, and Gauss point positions.
%
% REMARKS:
%   - This function is intended for use in problems where strain is assumed
%     small enough to justify the use of linearized kinematics.
%   - Unlike `BstLargeStrains`, the deformation gradient F is never computed.
%   - The output `IDENTITY_F` is left empty since F = I + J is not needed.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC - CIMNE, Barcelona
%   Date: 22-Nov-2020
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

% See  
% JAHO, 22-Nov-2020
if nargin ==0 
    load('tmp.mat')
end
% 
DATALOC.DEFORMATION_GRADIENT  = 0 ; 
MESH = DefaultField(MESH,'posgp_given',[]) ; 
MESH = DefaultField(MESH,'weights_given',[]) ; 

DATALOC.posgp_given = MESH.posgp_given ; 
DATALOC.weights_given = MESH.weights_given ; 
DATALOC.StrainStressWith4Components = MESH.DATA.StrainStressWith4Components ; 


[ Belem_F, wSTs, wST, XeALL, posgp] = ComputeBelemALL(MESH.COOR,MESH.CN,MESH.TypeElement,nstrain,DATALOC) ;
% Computing the deformation gradient as a vector 
% F = J + I , where J = Bst_F*d
[nelem,nnodeE ]= size(MESH.CN) ; 
ngaus = length(wSTs)/nelem ;
disp('Assembly of Bst (stacked B-matrix)...')
[nnode,ndim ]= size(MESH.COOR) ;
% if ndim == 2
%     nstrainF = 4;
%     IDENTITY = [1;1;0;0] ;
% elseif ndim ==3 
%     nstrainF = 9 ; 
%     IDENTITY = [1;1;1;0;0;0;0;0;0] ; 
% end

Bst_F = AssemblyBGlobal(Belem_F,nstrain,nelem,nnodeE,ndim,ngaus,MESH.CN,nnode) ;


IDENTITY_F  =  [] ; 


