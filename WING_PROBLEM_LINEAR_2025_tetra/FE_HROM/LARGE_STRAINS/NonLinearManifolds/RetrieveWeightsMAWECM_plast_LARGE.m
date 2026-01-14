function  [OPERFE,q_LAT] = RetrieveWeightsMAWECM_plast_LARGE(VAR,OPERFE)
%--------------------------------------------------------------------------
% RetrieveWeightsMAWECM_plast_LARGE
% ---------------------------------
% PURPOSE
%   Retrieve (and install in OPERFE) the MAW-ECM integration weights w and
%   their latent derivative dw/dq for a **plasticity / large-strain** setting,
%   evaluated at the current latent abscissa q_LAT taken from VAR.DISP.
%
% WHAT IT DOES
%   1) Extract latent coordinate used for MAW-ECM:
%        q_LAT = VAR.DISP(OPERFE.wSTs_cluster.IndexDOFl_q)
%   2) Evaluate the (precompiled) regressor for weights and gradient:
%        [w, dw_q] = nameFunctionEvaluate(q_LAT, DATA_regress)
%      where:
%        - w    : vector of ECM weights at the active latent point,
%        - dw_q : ∂w/∂q_LAT (same length as w).
%   3) Install results into OPERFE:
%        OPERFE.wSTs = w
%        OPERFE.dw_q = [0 ... 0, dw_q]   (pad to reduced space, placing dw_q
%                                          at column IndexDOFl_q)
%
% INPUTS
%   VAR    : State container; uses VAR.DISP as reduced coordinates.
%   OPERFE : Must contain wSTs_cluster struct with fields:
%            • IndexDOFl_q      : index of the latent DOF that parametrizes w
%            • DATA_regress     : {nameFunctionEvaluate, ...} descriptor for
%                                  MAW-ECM weights and their derivative.
%
% OUTPUTS
%   OPERFE : Updated with
%            • wSTs : current ECM weights at q_LAT
%            • dw_q : gradient of weights w.r.t. the latent coordinate,
%                     padded to full reduced-space width with zeros except at
%                     column IndexDOFl_q.
%   q_LAT  : The scalar latent abscissa used for evaluation.
%
% NOTES
%   • This helper assumes a **single controlling latent DOF** for the weight
%     law (the one pointed by IndexDOFl_q). For multi-latent control, extend
%     DATA_regress and generalize the padding to map each ∂w/∂q_j into the
%     corresponding column(s).
%   • w and dw_q dimensions must be consistent with the number of ECM points.
%   • Downstream assemblies (e.g., GetJacobianMatrix_MHROM) may use OPERFE.dw_q
%     when the **weight-term sensitivity** is enabled (MAW-ECM).
%
% DEPENDENCIES
%   • OPERFE.wSTs_cluster.DATA_regress.nameFunctionEvaluate(q, DATA_regress)
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

q_LAT = VAR.DISP(OPERFE.wSTs_cluster.IndexDOFl_q);
[w,dw_q]   = feval(OPERFE.wSTs_cluster.DATA_regress.nameFunctionEvaluate,q_LAT,OPERFE.wSTs_cluster.DATA_regress) ;
OPERFE.wSTs = w;
OPERFE.dw_q = [zeros(size(dw_q,1),OPERFE.wSTs_cluster.IndexDOFl_q-1),dw_q];