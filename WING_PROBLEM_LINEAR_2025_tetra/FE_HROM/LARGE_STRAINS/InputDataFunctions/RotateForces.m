function  Fbody =    RotateForces(Fbody,RotREFP,DATA,MESH)
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Rotate and compress a time history of nodal body forces to a (possibly
%   time-dependent) reference frame, then compute a low-rank basis via a
%   randomized SVD (RSVD). The routine:
%     1) applies a per-step rotation to the nodal force snapshots,
#     2) builds clustered snapshot matrices to control memory usage,
%     3) computes an RSVD over the rotated data,
%     4) updates Fbody.U (basis) and Fbody.a (reduced coordinates).
%
% INPUTS:
%   Fbody    : structure with snapshot factorization of forces.
%              • Fbody.U ∈ R^(ndof × r0) — current (or initial) basis.
%              • Fbody.a ∈ R^(r0 × Nsteps) — coefficients s.t. F ≈ U*a.
%              (If U,a come from a prior stage, they represent the forces
%               before rotation; this routine reprojects them after rotation.)
%   RotREFP  : 1×Nsteps cell array with rotation matrices for each time step.
%              Each RotREFP{t} ∈ R^(ndim × ndim) rotates vectors from the
%              global frame to the chosen rotating frame at step t.
%   DATA     : structure; uses DATA.MESH.ndof for memory estimation.
%   MESH     : structure with geometry.
%              • MESH.COOR ∈ R^(nnode × ndim) → ndof = nnode*ndim.
%
% OUTPUTS:
%   Fbody    : same structure with updated low-rank representation after
%              rotation:
%              • Fbody.U ∈ R^(ndof × r) — new basis from RSVD on rotated data.
%              • Fbody.a ∈ R^(r × Nsteps) — new coefficients (V'·S).
%              If RSVD returns empty U, falls back to zeros of compatible size.
%
% ALGORITHM OVERVIEW:
%   1) Dimensions and memory budget:
%      - ndof = nnode*ndim; Nsteps = size(Fbody.a,2).
%      - Estimate data size nsizeNOD = DATA.MESH.ndof*Nsteps*8e-6 (MB),
%        and split the time steps into 'nclusters' chunks so that each
%        cluster remains below LIMIT (MB).
%
%   2) Rotation and clustering:
%      - For each cluster c:
%          * Reconstruct forces at each step in the cluster:
%                Floc = Fbody.U * Fbody.a(:,t)  → (ndof×1)
%            reshape to nodal blocks (ndim × nnode), rotate with RotREFP{t}':
%                FlocROT = RotREFP{t}' * reshape(Floc,ndim,[])
%            and stack back into a column vector. Collect as columns of
%            RotatedMatrix_loc ∈ R^(ndof × nsteps_cluster).
%          * Store RotatedMatrix_loc in a cell array RotatedMatrix{c}.
%
%   3) Low-rank compression:
%      - Define EPSILON = 1e-10 (per cluster tolerance proxy) and call
%        RSVDqp(RotatedMatrix, EPSILON), which accepts a cell array of blocks
%        and returns (U,S,V) as if performing an SVD on the concatenation.
%      - If U is nonempty:
%            Fbody.U ← U,
%            Fbody.a ← (V' .* S)      % implemented via bsxfun(@times,V',S)
%        else:
%            Fbody.U ← zeros(ndof,1), Fbody.a ← zeros(1,Nsteps).
%
% IMPLEMENTATION NOTES:
%   - Memory control: clustering avoids forming a massive ndof×Nsteps matrix.
%   - Rotations: using RotREFP{t}' assumes column vectors represent components
%     in the global frame; transpose maps to the rotating frame.
%   - RSVDqp: expected to implement a numerically stable randomized SVD
%     over block inputs; S may be returned as a vector (singular values).
%
% ASSUMPTIONS / SANITY:
%   - RotREFP{t} must be orthogonal (rotation) of size ndim×ndim for all t.
%   - Fbody.U and Fbody.a must be conformable (size(U,1)=ndof, size(a,1)=size(U,2)).
%   - DATA.MESH.ndof equals ndof (used only for the memory estimate).
%   - No explicit checks on det(RotREFP{t}); if frames are not proper rotations,
%     magnitudes may be distorted.
%
% NUMERICAL CONSIDERATIONS:
%   - EPSILON is a small tolerance driving RSVD rank selection; adapt as needed.
%   - If forces are near-zero after rotation, RSVD may return empty U; code
%     handles this gracefully by zeroing U and a.
%
% REFERENCES:
%   - Internal note: “07_RotatingFrameStatic.pdf” (see code comment).
% -------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
end
nnode = size(MESH.COOR,1) ;
ndim = size(MESH.COOR,2) ;
ndof = nnode*ndim ;
LIMIT = 50 ; % ;Mb
% What is the size of nodal variables%
Nsteps = size(Fbody.a,2) ; 
nsizeNOD = DATA.MESH.ndof*Nsteps*8*1e-6 ;
nclusters = ceil(nsizeNOD/LIMIT);
NSTEPS_CLUSTER = cell(1,nclusters) ;

% We have to divide  1:Nsteps into nclusters
FREQ  = ceil(Nsteps/nclusters) ;
iini = 1;
%ifin = FREQ ;
for i=1:nclusters
    ifin = iini + FREQ-1 ;
    ifin = min(ifin,Nsteps) ;
    NSTEPS_CLUSTER{i} = iini:ifin ;
    iini = ifin +1 ;
end
RotatedMatrix = cell(1,nclusters) ;
for icluster = 1:nclusters
    nstepsLOC = length(NSTEPS_CLUSTER{icluster}) ;
    RotatedMatrix_loc = zeros(ndof,nstepsLOC);
    for istepLOC = 1:nstepsLOC
        istep = NSTEPS_CLUSTER{icluster}(istepLOC) ;
        Floc = Fbody.U*Fbody.a(:,istep) ;
        Floc = reshape(Floc,ndim,[]) ;
        FlocROT = RotREFP{istep}'*Floc ;
        RotatedMatrix_loc(:,istepLOC) = FlocROT(:) ;
    end
    RotatedMatrix{icluster} = RotatedMatrix_loc ;
end

EPSILON  = 1e-10*ones(1,nclusters) ;

[U,S,V] = RSVDqp(RotatedMatrix,EPSILON) ;

if  ~isempty(U)
Fbody.U = U ;
a = bsxfun(@times,V',S) ;
Fbody.a = a;
else
    Fbody.U = zeros(ndof,1) ; 
    Fbody.a = zeros(1,Nsteps) ; 
end

 
