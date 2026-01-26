function d = ObtainDisplacementFromSnapshotsProj(BasisU,qL_extended,SNAP_cluster,DOFr,OTHER_output,DOFl)
%==========================================================================
% ObtainDisplacementFromSnapshotsProj
%
% PURPOSE
%   Reconstruct the FULL displacement field d from:
%     • a reduced (projected) approximation of the FREE DOFs, dL ≈ BasisU*qL_extended,
%     • the corresponding reference (SVD) reconstruction of the PRESCRIBED DOFs, dR_exact,
%   and (if present) an **affine Dirichlet coupling** d = A*dL + [0; G*d_m] that links
%   slave/prescribed boundary DOFs to a subset of master DOFs. The routine enforces
%   essential BCs in the weak-form sense, consistent with the slides’ treatment of
%   trial/test spaces, block partitioning, and vectorized assembly.
%
% THEORETICAL CONTEXT (affine/essential BCs)
%   • In the weak form, essential (Dirichlet) conditions are enforced on the trial space
%     and carried to the discrete level by splitting dofs into free (l) and prescribed (r)
%     blocks and/or by **eliminating** slave DOFs via an affine map. :contentReference[oaicite:1]{index=1}
%   • Classic block form (linear reference) is
%         [Krr Krl; Klr Kll][dr; dl] = [Fr; Fl]  ⇒  Klldl = Fl − Klrdr,
%     highlighting that only dl are unknowns while dr are known/prescribed. :contentReference[oaicite:2]{index=2}
%   • In manifold HROM with moving (non-homogeneous) Dirichlet data, boundary DOFs are
%     handled by an **affine constraint**
%         d = A*dL  +  [0; G*d_m],
%     where A injects free dofs into the full vector and G maps a subset of master dofs
%     (DOFm) onto the prescribed/slave set DOFr. This matches the slides’ “essential vs
%     natural” BCs and the vectorized assembly viewpoint (B, L, N operators).
%
% WHAT THIS ROUTINE DOES
%   1) Rebuild the local/free approximation:  dL = BasisU*qL_extended.
%   2) Recover an **exact** snapshot-based reconstruction of the prescribed part:
%        dR_exact = U_r * S * V^T   on DOFr
%      (taken from the SVD factors stored in SNAP_cluster), ensuring that the final d
%      is consistent with the offline snapshots at boundary nodes. This is coherent with
%      the weak-form imposition of Dirichlet data (prescribed “r” dofs). :contentReference[oaicite:4]{index=4}
%   3) If NO affine BC operator A is provided:
%        – Assemble d by placing dL at DOFl and dR_exact at DOFr.
%   4) If an affine BC operator A **is** provided (OTHER_output.DISP_CONDITIONS.A):
%        – Start from d = A*dL (inject free approx into the full vector).
%        – Because the affine law is d = A*dL + CONSTANT_TERM, add the missing
%          **boundary constant** inferred from the exact snapshot:
%               • Build dEXACT by inserting the exact local part dL_exact at DOFl
%                 from the snapshot SVD.
%               • Using the master set DOFm and the coupling matrix G (in
%                 OTHER_output.DISP_CONDITIONS.G), compute
%                     ConstantTERM_R = dR_exact − G * dEXACT(DOFm,:)
%                 and add it to d on DOFr.
%      This guarantees that the reconstructed full field d matches the snapshot’s
%      prescribed boundary motion, while the interior respects the reduced
%      approximation—consistent with essential BC enforcement in the slides. :contentReference[oaicite:5]{index=5}
%
% WHY THIS IS CONSISTENT WITH THE SLIDES
%   • Essential BCs: imposed in the trial space (functions satisfy u(0)=g in the 1D
%     model), hence prescribed dofs are not solved for but reconstructed/inserted. :contentReference[oaicite:6]{index=6}
%   • Block logic & elimination: free vs prescribed sets mirror the Kll/Klr structure. :contentReference[oaicite:7]{index=7}
%   • Vectorized operators: the full d feeds B, L, N to compute strains/forces
%     consistently at Gauss points in subsequent steps. :contentReference[oaicite:8]{index=8}
%
% INPUTS
%   BasisU        : Reduced basis spanning the free/local displacements (columns).
%   qL_extended   : Reduced coordinates after encoder application (may include nonlinear
%                   complement from τ(q)); size = n_modes × n_snapshots.
%   SNAP_cluster  : Structure holding SVD factors for displacement snapshots:
%                   .DISP.U, .DISP.S, .DISP.V (partitioned by DOFl/DOFr).
%   DOFr          : Global indices of prescribed/slave (Dirichlet) DOFs.
%   OTHER_output  : May contain DISP_CONDITIONS with fields:
%                      .A    (injector from free to full),
%                      .G    (slave–master coupling on boundary),
%                      .DOFm (master DOF indices participating in G).
%   DOFl          : Global indices of free DOFs.
%
% OUTPUT
%   d             : Full-order displacement snapshot(s), size (ndof × n_snapshots),
%                   consistent with affine Dirichlet constraints when present.
%
% NOTES / PITFALLS
%   • If A and G are inconsistent with the partition (DOFl/DOFr/DOFm), boundary
%     corrections may be misapplied. Ensure DOFm ⊂ DOFl ∪ boundary-masters.
%   • Using dR_exact to infer the constant term makes reconstructed d exactly match
%     the boundary motion of the training snapshot—this is desirable for stress and
%     reaction consistency in later offline/online steps.
% JAHO, 8-Oct-2025, Balmes 185, BARcelona
% Comments by ChatGPT5
%==========================================================================
if nargin == 0
    load('tmp1.mat')
end


dL = BasisU*qL_extended; % Snapshot displacements (local)
dR_exact = bsxfun(@times,SNAP_cluster.DISP.U(DOFr,:)',SNAP_cluster.DISP.S)' ;
dR_exact = dR_exact*SNAP_cluster.DISP.V' ;
ndof = size(dL,1)+size(dR_exact,1) ;
if isempty(OTHER_output.DISP_CONDITIONS.A)
    
    d = zeros(ndof,size(dL,2)) ;
    
    
    d(DOFl,:)  = dL ;
    if size(dL,2) == size(dR_exact,2)
    d(DOFr,:)  = dR_exact ;
    else
         d(DOFr,:)  = dR_exact(:,2:end) ; 
    end
else
    % AFfine boundary conditions
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/13_HOMOG.mlx
    d = OTHER_output.DISP_CONDITIONS.A*dL ;
    % In this case d = A*dL_approx + CONSTANT_TERM
    % However, we do not have at our disposal CONSTANT_TERM
    % Let us infer it from the "exact" displacement
    dL_exact = bsxfun(@times,SNAP_cluster.DISP.U(DOFl,:)',SNAP_cluster.DISP.S)' ;
    dL_exact = dL_exact*SNAP_cluster.DISP.V' ;
    dEXACT = zeros(ndof,size(dL,2)) ;
    dEXACT(DOFl,:) = dL_exact ;
    DOFm = OTHER_output.DISP_CONDITIONS.DOFm ;
    ConstantTERM_R = dR_exact-OTHER_output.DISP_CONDITIONS.G*dEXACT(DOFm,:);
    d(DOFr,:) = d(DOFr,:) + ConstantTERM_R ;
end