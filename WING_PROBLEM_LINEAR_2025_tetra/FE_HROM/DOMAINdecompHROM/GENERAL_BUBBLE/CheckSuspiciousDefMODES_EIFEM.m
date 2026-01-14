function PhiDEFbs_INV = CheckSuspiciousDefMODES_EIFEM(DATAcommon,Sbs,DATAoffline,PhiDEFbs,MESH)
%--------------------------------------------------------------------------
% function PhiDEFbs_INV = CheckSuspiciousDefMODES_EIFEM( ...
%                          DATAcommon, Sbs, DATAoffline, PhiDEFbs, MESH)
%
% PURPOSE:
%   This function analyzes the spectrum of singular values associated with 
%   deformation modes (Φ_def) to detect *invariant* modes, i.e., modes that
%   are preserved across variations in geometry or loading conditions.
%
%   It either automatically selects these invariant deformation modes based 
%   on singular value decay (heuristic thresholding) or loads them from an 
%   externally provided dataset (e.g., if defined over a different mesh).
%
% INPUTS:
% -------
%   - DATAcommon: structure with global/common settings. If it contains
%       · DATAcommon.PhiDEFbs_GIVEN_OTHER_MESH → uses external invariant modes.
%
%   - Sbs: vector of singular values associated with deformation modes
%
%   - DATAoffline: offline settings, including:
%       · TOLSVD_INVARIANT_BASIC_MODES
%       · NumberOfInvariantModes_SVD
%       · THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES
%
%   - PhiDEFbs: [nDOFs × nModes] matrix of all deformation modes (Φ_def)
%
%   - MESH: current mesh structure (used to verify conformity with stored modes)
%
% OUTPUT:
% -------
%   - PhiDEFbs_INV: [nDOFs × nInvModes] matrix of selected *invariant* 
%     deformation modes.
%
% METHOD:
% -------
% CASE 1: No predefined invariant basis (automatic selection):
%   → The function computes the ratios of consecutive singular values:
%
%         RATIO_SV(i) = Sbs(i) / Sbs(i+1)
%
%   → If a sharp drop is detected (ratio > threshold), the modes before
%     the drop are selected as invariant.
%
%   → Otherwise, it selects the number of invariant modes set in
%     DATAoffline.NumberOfInvariantModes_SVD.
%
% CASE 2: External invariant basis provided (e.g., from different mesh):
%   → Performs mesh centering and KNN search to match nodes
%   → Transfers basis modes to current mesh using a degree-of-freedom mapping
%
% THEORY:
% -------
%   This procedure supports the Empirical Invariant Finite Element Method (EIFEM),
%   where a portion of the reduced basis is *universal* (e.g., from unit cell training)
%   and the rest is enriched online.
%
%   The logic is tied to the spectral decay of deformation snapshots and assumes
%   that invariant modes dominate the first part of the spectrum.
%
%   If the singular value decay is mild (no significant drop), the algorithm defaults
%   to user-specified truncation.
%
% EXAMPLE:
% -------
%   [Φ_inv] = CheckSuspiciousDefMODES_EIFEM(DATAcommon, Sbs, DATAoffline, Φ, MESH);
%
% NOTES:
% ------
%   - Mesh conformity is required if modes are transferred from another mesh.
%   - Uses `knnsearch` and `small2large` for spatial DOF mapping.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE, May 2024
%--------------------------------------------------------------------------


if isempty(DATAcommon.PhiDEFbs_GIVEN_OTHER_MESH)
    
    disp(['Singular Values Deformational Modes (in the search for invariant modes)'])
    
    Sbs/Sbs(1)
    
    RATIO_SV = Sbs(1:end-1)./Sbs(2:end) ;
    
    if DATAoffline.TOLSVD_INVARIANT_BASIC_MODES >0
        THinv = 1e20 ;
    else
        THinv = 100 ;
    end
    
    DATAoffline = DefaultField(DATAoffline,'THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES',THinv) ;
    
    
    III =   find(RATIO_SV>DATAoffline.THRESHOLD_RATIO_CONSECUTIVE_SINGLE_VALUES_INVARIANT_DEF_MODES) ;
    
    if isempty(III)
        nINVmodes = length(Sbs) ;
    else
        nINVmodes = III(1) ;
    end
    disp('***********************************************************************************************************')
    disp('Searching for invariant modes')
    disp(['Total number of  deformational modes (for the first set of training tests) = ',num2str(length(Sbs))]) ;
    disp(['Number of  deformational modes that may be deemed INVARIANT = ',num2str(nINVmodes)]) ;
    if ~isempty(III)
        disp(['Jump in singular values (ratio consecutive singular values) =',num2str(RATIO_SV(nINVmodes))])
        disp('***********************************************************************************************************')
    else
        disp(['tolerance SVD = ',num2str(DATAoffline.TOLSVD_INVARIANT_BASIC_MODES)])
    end
    
    
    nINVmodes_select = min(nINVmodes,DATAoffline.NumberOfInvariantModes_SVD) ; 
    
    if  nINVmodes == DATAoffline.NumberOfInvariantModes_SVD
    disp(['Final number of modes selected by the user (variable DATAoffline.NumberOfInvariantModes SVD) = ',num2str(DATAoffline.NumberOfInvariantModes_SVD)])
     
    end
    
    
    PhiDEFbs_INV = PhiDEFbs(:,1:nINVmodes_select) ;
    
    
else
    disp('***************************************************************++')
    disp(['Invariant deformational modes given by the user'])
    disp('***************************************************************++')
    COORold = DATAcommon.PhiDEFbs_GIVEN_OTHER_MESH.COOR;
    if size(COORold,1) ~= size(MESH.COOR,1)
        error('The number of nodes (and their locations ) should be identical')
    end
    COORnew = MESH.COOR;
    COORnew_cent = sum(COORnew,1)/size(COORnew,1) ;
    COORold_cent = sum(COORold,1)/size(COORold,1) ;
    
    for idim = 1:size(COORold,2)
        COORnew(:,idim) = COORnew(:,idim)-COORnew_cent(idim) ;
        COORold(:,idim) = COORold(:,idim)-COORold_cent(idim) ;
    end
    
    disp(['Checking whether the meshes are conforming, and perform the spatial search'])
    
    [IdxNEW,DDD]  = knnsearch(COORnew,COORold) ;
    
    if max(DDD) >1e-10
        error('Not conforming meshes')
    end
    
    IdxDOFS = small2large(IdxNEW,size(COORnew,2)) ;
    
    PhiDEFbs_INV = DATAcommon.PhiDEFbs_GIVEN_OTHER_MESH.PhiDEFbs(IdxDOFS,:);
    
    
end