function    [OPERFE,MESH,DATA] = FEoperators_1dom(DATA,MESH,MATPRO)
%--------------------------------------------------------------------------
% FUNCTION: FEoperators_1dom
%
% PURPOSE:
%   Initializes and computes finite element (FE) operators and related mesh
%   quantities for a given domain. This includes:
%     - Normals and tangents to the boundary surfaces (for Neumann BCs or fluxes),
%     - Gauss point locations and weights,
%     - B-matrix (strain-displacement matrix),
%     - Mass matrix (if applicable),
%     - Metadata for stress/strain evaluations at Gauss points.
%
% USAGE:
%   [OPERFE, MESH, DATA] = FEoperators_1dom(DATA, MESH, MATPRO)
%
% INPUT:
%   - DATA   : Structure containing general simulation and discretization parameters.
%   - MESH   : Structure containing node coordinates, connectivities, and optionally boundary entities.
%   - MATPRO : Material property structure (may be used internally by `GeometricMatricesFunMULT`).
%
% OUTPUT:
%   - OPERFE : Structure containing precomputed FE operators such as:
%                * Bst  : strain-displacement matrix at all Gauss points
%                * wSTs : weights of Gauss integration
%                * posgp: positions of Gauss points
%                * M    : mass matrix (if computed)
%   - MESH   : Updated mesh structure with tangent/normal information and
%             assigned Gauss point data.
%   - DATA   : Updated structure with derived quantities such as:
%                * ngausT           : total Gauss points in domain
%                * ndofSTRESS       : degrees of freedom for stress field
%                * ndim             : number of spatial dimensions
%                * ndof             : total DOFs (nnode * ndim)
%
% ALGORITHM OVERVIEW:
%   1. Compute normal and tangent vectors on mesh boundaries (via NormalsAndTangents).
%   2. Assign or propagate Gauss point positions and weights (from `DATA` if given).
%   3. Compute core FE matrices and operators using `GeometricMatricesFunMULT`.
%   4. Populate internal fields of `DATA.MESH` with key quantities for further use
%      in assembly and integration.
%
% REFERENCES:
%   - Developed as part of a reduced-order modeling framework (EIFEM).
%   - Inspired by: 
%     /FE_HROM/LARGE_STRAINS/InputDataFunctions/PreProcessInputDataDyn1.m
%   - See test case: 
%     /TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
%
% AUTHOR:
%   Joaquín A. Hernández Ortega
%   UPC-CIMNE
%   Date: 9-Feb-2023
%   Comments by ChatGPT4, 13-May-2025
%--------------------------------------------------------------------------

% Function patterned after /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/FE_HROM/LARGE_STRAINS/InputDataFunctions/PreProcessInputDataDyn1.m
% It computes some FE operators and variables, such as the mass matrix, and
% the B-matrix (strain-displacement) operator 
% JAHO, 9-Feb-2023
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx

%
 if nargin == 0
    load('tmp.mat')
end
% ----------------------------
%%%%%%%%%%%%%%%%%
% ---------------
% NORMALS AND TANGENT VECTORS TO EACH BOUNDARY SURFACE
% ----------------------------------------------
MESH = NormalsAndTangents(MESH) ; 
 
[nnode,ndim ]= size(MESH.COOR) ;% Number of nodes
[nelem,nnodeE ]= size(MESH.CN) ; % Number of elements

DATA = DefaultField(DATA,'posgp_given',[]) ; 
DATA = DefaultField(DATA,'weights_given',[]) ; 
MESH.posgp_given =DATA.posgp_given ; 
MESH.weights_given =DATA.weights_given ; 


% 4. FE OPERATORS,MATRICES
% ---------------------------
MESH.DATA = DATA; 
% FINITE ELEMENT OPERATORS, AS WELL AS OTHER PROPERTIES OF THE DOMAIN (MASS, CENTROID...)
[OPERFE,MESH] =  GeometricMatricesFunMULT(MESH,MATPRO) ; 

 
DATA.MESH.posgp  =OPERFE.posgp ;
DATA.MESH.ngaus  =OPERFE.ngaus_STRESS ;

 
DATA.MESH.ngausT =  length(OPERFE.wSTs) ; % Total number of Gauss points
DATA.MESH.ndofSTRESS =   DATA.MESH.ngausT*DATA.MESH.nstrain; % Total number of Gauss points
DATA.MESH.ngaus_STRESS = OPERFE.ngaus_STRESS;
DATA.MESH.ndof = nnode*ndim;
DATA.MESH.ndim = ndim;

 
% 