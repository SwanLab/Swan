function   Fpnt =  GaussianTraction(nnodeALL,ndim,ISLOCAL,FORCE_INFO,NstT_W,MESH,iface) ; 
%==========================================================================
% function Fpnt = GaussianTraction(nnodeALL, ndim, ISLOCAL, FORCE_INFO, ...
%                                   NstT_W, MESH, iface)
%--------------------------------------------------------------------------
% Purpose
% -------
% Computes the nodal force vector associated with a Gaussian-type
% Neumann traction applied over a given boundary face of the mesh.
%
% The load represents a spatially localized traction (e.g. vertical
% pressure patch) whose magnitude follows a Gaussian distribution
% along one local coordinate direction of the surface.
%
% Input arguments
% ----------------
% nnodeALL : (integer) total number of nodes in the full mesh.
%
% ndim     : (integer) spatial dimension (2 or 3).
%
% ISLOCAL  : (logical flag) indicates whether the traction is defined in
%             local or global coordinates. (Not used explicitly here but
%             kept for interface compatibility.)
%
% FORCE_INFO : structure containing Gaussian traction parameters:
%     .AMPLITUDE        -> Vector [Tx, Ty, Tz] with the traction
%                          components (peak values) [N/m²].
%     .CENTER           -> Center of the Gaussian, expressed in the
%                          selected coordinate direction (scalar).
%     .SIGMA_WIDTH      -> Gaussian width parameter σ (in same units
%                          as coordinates, e.g. meters).
%     .COORDINATE_LOCAL -> Index of the coordinate direction along which
%                          the Gaussian varies (1=x, 2=y, 3=z).
%
% NstT_W : (unused placeholder, kept for interface compatibility with
%           other traction definitions or time weighting schemes).
%
% MESH : structure containing geometric and topological data:
%     .NODES_FACES{iface}         -> Node indices of face 'iface'.
%     .COOR                       -> Nodal coordinates matrix [nnode x ndim].
%     .PROPERTIES_FACES{iface}.GeometricMassMatrix
%                                 -> Local geometric (surface) mass matrix
%                                    for the given face, used to integrate
%                                    the traction over the element surface.
%
% iface : (integer) index of the face on which the load is applied.
%
% Output argument
% ----------------
% Fpnt : Global force vector [nnodeALL*ndim x 1], containing the nodal
%        equivalent forces due to the Gaussian traction on face 'iface'.
%
%--------------------------------------------------------------------------
% Computational steps
% -------------------
% 1. Extract the coordinates of the nodes on the selected face.
% 2. Compute the Gaussian exponent along the chosen coordinate direction:
%        eEXP = exp( - (x - CENTER)^2 / (2*SIGMA_WIDTH^2) )
%
% 3. Evaluate the traction magnitude at each node as
%        t(x) = AMPLITUDE .* eEXP
%
% 4. Multiply by the face geometric mass matrix (Mst) to obtain the
%    consistent nodal force contributions:
%        F_face = Mst * t(x)
%
% 5. Scatter these local contributions into the global nodal force vector
%    Fpnt, using the appropriate degrees of freedom for each dimension.
%
%--------------------------------------------------------------------------
% Physical interpretation
% -----------------------
% The traction applied on the face has the spatial form:
%     t_i(x_local) = AMPLITUDE_i * exp[-(x_local - CENTER)^2 / (2σ²)]
%
% The resulting total (integrated) force per component is approximately
%     F_i ≈ AMPLITUDE_i * √(2π) * σ,
% provided that the face length is much larger than σ.
%
%--------------------------------------------------------------------------
% Notes
% -----
% • The current implementation applies the traction independently to
%   each Cartesian component defined in AMPLITUDE.
%
% • The geometric mass matrix Mst should already include the appropriate
%   Gauss integration weights and Jacobian determinants for the surface.
%
% • This function is typically called by Neumann assembly routines when
%   evaluating time-dependent or moving Gaussian load patterns.
%
%--------------------------------------------------------------------------
% Example
% -------
%  FORCE_INFO.AMPLITUDE = [0, -1e6];    % vertical downward 1 MPa
%  FORCE_INFO.SIGMA_WIDTH = 0.00125;    % 1.25 mm
%  FORCE_INFO.CENTER = 0.5;             % centered at x = 0.5 m
%  FORCE_INFO.COORDINATE_LOCAL = 1;     % varies along x_local
%
%  Fpnt = GaussianTraction(nnode, 2, true, FORCE_INFO, [], MESH, 3);
%
%--------------------------------------------------------------------------
% J.A. Hernández-Ortega (JAHO)
% Starbucks, Paseo de Gracia, Barcelona
% 8-Nov-2025 (Saturday)
%==========================================================================

if nargin == 0
  load('tmp.mat')  
end

%
nodesfA  = MESH.NODES_FACES{iface} ; 

COOR_FACE = MESH.COOR(nodesfA,FORCE_INFO.COORDINATE_LOCAL) ; % Coordinates of this face
DeltaX_2=  (COOR_FACE - FORCE_INFO.CENTER).^2 ;
Exponente = DeltaX_2/(2*FORCE_INFO.SIGMA_WIDTH^2) ; 
eEXP = exp(-Exponente) ; 

Mst = MESH.PROPERTIES_FACES{iface}.GeometricMassMatrix  ; 

Fpnt = zeros(size(MESH.COOR,1)*ndim,1) ;
ndim  =size(MESH.COOR,2) ; 
    DOFsGLO = small2large(nodesfA,ndim) ; 

tMAX = FORCE_INFO.AMPLITUDE(:) ; 
for idim  =1:ndim
    Fpnt_loc = Mst*tMAX(idim)*eEXP ;  
    DOFsGLO_i = DOFsGLO(idim:ndim:end) ; 
    Fpnt(DOFsGLO_i) =  Fpnt_loc ; 
end

% 
%  






 