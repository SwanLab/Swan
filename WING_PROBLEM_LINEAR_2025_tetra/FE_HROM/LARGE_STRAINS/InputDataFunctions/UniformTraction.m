function   Fpnt =  UniformTraction(nnodeALL,ndim,ISLOCAL,FORCE_INFO,NstT_W,MESH,iface) ; 
%==========================================================================
% function Fpnt = UniformTraction(nnodeALL, ndim, ISLOCAL, FORCE_INFO, ...
%                                 NstT_W, MESH, iface)
%--------------------------------------------------------------------------
% Purpose
% -------
% Assemble the global nodal force vector produced by a **uniform** Neumann
% traction acting on a *subset* of a boundary face. This routine mirrors
% GaussianTraction but replaces the Gaussian footprint by a **top-hat**
% (constant) footprint of finite length along a chosen surface coordinate.
%
% The active segment is centered at FORCE_INFO.CENTER and has length
% FORCE_INFO.SIGMA_WIDTH (interpreted here as the footprint length “a”).
%
%--------------------------------------------------------------------------
% Inputs
% ------
% nnodeALL : (int)  total number of mesh nodes (size(MESH.COOR,1)).
%
% ndim     : (int)  spatial dimension (2 or 3). Used to expand scalar
%                    nodal forces into vector DOFs (x,y,(z)).
%
% ISLOCAL  : (logical) kept for interface compatibility (not used here).
%
% FORCE_INFO : struct with uniform-traction parameters
%   .AMPLITUDE        -> [Tx Ty (Tz)] components of traction (N/m^2),
%                        per unit out-of-plane thickness in plane strain.
%   .CENTER           -> center of the uniform footprint along the chosen
%                        local coordinate (same units as coordinates).
%   .SIGMA_WIDTH      -> footprint length “a” (NOT a Gaussian σ here);
%                        the load acts on [CENTER - a/2, CENTER + a/2].
%   .COORDINATE_LOCAL -> index of coordinate used to define the segment
%                        along the face: 1=x_local, 2=y_local, (3=z_local).
%
% NstT_W  : (placeholder) time weighting / extra arguments (unused).
%
% MESH    : struct with boundary topology and geometry
%   .NODES_FACES{iface}                    -> node IDs on face 'iface'
%   .Indexes_faces_bnd_element{iface}      -> indices of boundary elements
%   .CNb                                   -> boundary-element connectivities
%   .COOR                                  -> nodal coordinates [nnode x ndim]
%   .TypeElementB                          -> boundary element type/tag
%
% iface   : (int) boundary face index where the traction is applied.
%
%--------------------------------------------------------------------------
% Output
% ------
% Fpnt : (nnodeALL*ndim x 1) global vector of equivalent nodal forces due
%        to the **uniform** traction over the selected segment of face 'iface'.
%
%--------------------------------------------------------------------------
% Method / Algorithm
% ------------------
% 1) Identify nodes belonging to face 'iface':
%       nodesfA = MESH.NODES_FACES{iface}
%
% 2) Select the **active subset** of face nodes whose chosen coordinate
%    lies in the interval:
%       x ∈ [CENTER - (SIGMA_WIDTH/2)*TOL,  CENTER + (SIGMA_WIDTH/2)*TOL]
%    where TOL = 1+1e-6 provides a tiny robustness margin.
%
% 3) Determine the boundary elements “under” the active node subset:
%       Ind_facesElementsBND = MESH.Indexes_faces_bnd_element{iface}
%       CONNECTb = MESH.CNb(Ind_facesElementsBND,:)
%       [CNbLOC, setBelem] = ElemBnd(CONNECTb, nodesfA_subset_glo)
%
% 4) Build a **geometric (surface) mass matrix** Mst for those boundary
%    elements (length/area integration) and form local uniform nodal loads:
%       [CentroidFA, AREA, Mst] = CentroidGeometricMassMatrixNEW(...)
%       LocalForceNodal = Mst * 1     // “1” = uniform traction shape
%
% 5) Scatter the scalar nodal loads into the global vector by components:
%       GlobalForceNodal(nodesfA_subset_glo) = LocalForceNodal
%       For each idim = 1:ndim:
%           Fpnt(idim:ndim:end) += tMAX(idim) * GlobalForceNodal
%    where tMAX = FORCE_INFO.AMPLITUDE(:).
%
%--------------------------------------------------------------------------
% Physical interpretation
% -----------------------
% • The traction is **constant** over the active segment and zero elsewhere:
%       t_i = AMPLITUDE_i,   for x ∈ [CENTER - a/2, CENTER + a/2]
%       t_i = 0,             otherwise
%   with a = SIGMA_WIDTH in this routine (top-hat width).
%
% • Resultant (per unit thickness) along component i is approximately:
%       F_i ≈ AMPLITUDE_i * (length of active segment on the face).
%   The geometric mass matrix Mst ensures consistent nodal distribution
%   according to element metrics and quadrature.
%
%--------------------------------------------------------------------------
% Conventions & Units
% -------------------
% • AMPLITUDE is in N/m^2 (line traction per unit thickness in plane strain).
% • Coordinates, CENTER, SIGMA_WIDTH in meters.
% • The routine acts only on the specified face 'iface'; ensure the face
%   parameterization uses COORDINATE_LOCAL monotonically along the segment.
%
%--------------------------------------------------------------------------
% Notes / Caveats
% ---------------
% • This function interprets SIGMA_WIDTH as a **uniform footprint length**,
%   unlike GaussianTraction where it is a standard deviation σ.
%
% • The tolerance TOL = 1+1e-6 slightly expands the selection interval to
%   avoid losing nodes on floating-point boundaries; adjust if needed.
%
% • The helper routines ElemBnd and CentroidGeometricMassMatrixNEW must be
%   consistent with MESH.TypeElementB and the face connectivity ordering.
%
% • For time-dependent or moving loads, update CENTER (and possibly the
%   active face index) across time steps and call this routine per step.
%
%--------------------------------------------------------------------------
% Example
% -------
% FORCE_INFO.AMPLITUDE       = [0; -1e6];   % 1 MPa downward
% FORCE_INFO.CENTER          = 0.50;        % center at x = 0.5 m
% FORCE_INFO.SIGMA_WIDTH     = 0.005;       % active footprint a = 5 mm
% FORCE_INFO.COORDINATE_LOCAL= 1;           % vary along local x
% Fpnt = UniformTraction(nnode, 2, true, FORCE_INFO, [], MESH, 3);
%
%--------------------------------------------------------------------------
% J.A. Hernández-Ortega (JAHO), 9-Nov-2025, Barcelona
%==========================================================================


if nargin == 0
  load('tmp.mat')  
end

%
nodesfA  = MESH.NODES_FACES{iface} ; % These are the nodes of the surface that is being transversed by the uniform load

COOR_FACE = MESH.COOR(nodesfA,FORCE_INFO.COORDINATE_LOCAL) ; % Coordinate along it moves, as many as nodes
TOL = 1+1e-6; 
xINTERVAL =   [FORCE_INFO.CENTER-(FORCE_INFO.SIGMA_WIDTH/2)*TOL,FORCE_INFO.CENTER+(FORCE_INFO.SIGMA_WIDTH/2)*TOL] ; 
nodesfA_subset_loc = find(COOR_FACE<= xINTERVAL(2) & COOR_FACE>= xINTERVAL(1)) ;  % Subset of nodes of the elements under consideration


nodesfA_subset_glo = nodesfA(nodesfA_subset_loc) ; 

% But which are the actual elements ?  What are their connectivities ? 
Ind_facesElementsBND = MESH.Indexes_faces_bnd_element{iface} ; % Indexes boundary elements
% Local connectivities, boundary elements
CONNECTb = MESH.CNb(Ind_facesElementsBND,:) ; 
[CNbLOC, setBelem]= ElemBnd(CONNECTb,nodesfA_subset_glo) ; 
 [CentroidFA,AREA,Mst] =CentroidGeometricMassMatrixNEW(MESH.COOR,[],CNbLOC,MESH.TypeElementB) ;
% Local forces, uniform 
LocalForceNodal  = Mst*ones(length(nodesfA_subset_glo),1) ; 
GlobalForceNodal = zeros(size(MESH.COOR,1),1) ; 
GlobalForceNodal(nodesfA_subset_glo) = LocalForceNodal ; 

 
 
Fpnt = zeros(size(MESH.COOR,1)*ndim,1) ;
ndim  =size(MESH.COOR,2) ; 

tMAX = FORCE_INFO.AMPLITUDE(:) ; 
for idim  =1:ndim
     
    Fpnt(idim:ndim:end) =  GlobalForceNodal*tMAX(idim) ; 
end

% 
%  






 