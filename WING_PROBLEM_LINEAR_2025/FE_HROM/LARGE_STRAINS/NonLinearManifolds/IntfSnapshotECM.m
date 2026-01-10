function A_fint = IntfSnapshotECM(qL_extended,nREDcoor,DATA_evaluateTAU_and_DER,BstRED_l,SNAPstressSTWOproj_LOC_INELAST,DATA)
%--------------------------------------------------------------------------
% IntfSnapshotECM
%
% PURPOSE:
%   Assemble the nonlinear internal force snapshots needed for the ECM 
%   hyperreduction procedure. The function evaluates, at each generalized 
%   coordinate qL, the reduced internal force contribution associated with 
%   nonlinear stresses rather than the raw stresses.
%
% INPUTS:
%   qL_extended : [nq × nt] matrix of reduced generalized coordinates 
%                 for all snapshots (nq = #reduced DOFs, nt = #snapshots).
%
%   nREDcoor    : Number of reduced generalized coordinates (dimension of q).
%
%   DATA_evaluateTAU_and_DER : Structure with fields specifying the decoder
%                 τ(q) and its derivatives, including the handle to the 
%                 evaluation function (nameFunctionEvaluate).
%
%   BstRED_l    : Reduced operator Bst (linearized strain-displacement operator
%                 projected onto reduced basis).
%
%   SNAPstressSTWOproj_LOC_INELAST : [nsnap × nt] matrix of projected 
%                 nonlinear stress snapshots (e.g. PK1 or PK2 stresses).
%
%   DATA        : Structure containing additional problem data and options 
%                 (e.g. quadrature weights, constitutive flags).
%
% OUTPUT:
%   A_fint      : Cell array of size {1 × nt}, where each entry contains
%                 the reduced internal force operator (integrand form) 
%                 at the corresponding snapshot.
%
% CONTEXT:
%   Used in nonlinear manifold HROM/ECM workflows where internal forces are
%   reconstructed from stress snapshots and τ′(q). The integrand A_fint 
%   forms the building block for computing reduced internal forces and 
%   tangent matrices.
%
% JAHO — August 2025
%--------------------------------------------------------------------------


% Now we use nonlinear stresses rather than the actual stresses !
A_fint = cell(1,size(qL_extended,2));                 % one per time step
% nAfint = zeros(1,size(qL_extended,2)) ;
for itimeLOC = 1:size(qL_extended,2)
    q = qL_extended(1:nREDcoor,itimeLOC);             % generalized coords
    % τ′(q): needed to map Bst*Φ*τ′(q) into current tangent
    [~,tauNONder_q,~] = feval(DATA_evaluateTAU_and_DER.nameFunctionEvaluate, ...
        q, DATA_evaluateTAU_and_DER);
    BstRED_l_q = BstRED_l * tauNONder_q;              % reduced B at current q
    %   Pk1_stress = SNAPstressPonePROJ_LOC{iloc}(:,itimeLOC); % PK1 vector at t
    Stress_nonlinear = SNAPstressSTWOproj_LOC_INELAST(:,itimeLOC) ;
    % Build A_fint column (internal force integrand mapped by PK1)
    A_fint{itimeLOC} = BasisF_from_BasisStress_PK1(BstRED_l_q,Stress_nonlinear,DATA);
    
    %  nAfint(itimeLOC) = norm(A_fint{itimeLOC},'fro');
    
end