function  [OPERFE,qNONrel] = RetrieveWeightsMAWECM_DAMAGE(VAR,OPERFE)
%--------------------------------------------------------------------------
% RetrieveWeightsMAWECM_DAMAGE
% ---------------------------------
% PURPOSE
%   Compute MAW-ECM integration weights for a DAMAGE model in which the
%   weight law depends on a **ratio of two latent coordinates**:
%       qLIN  = linear driving latent DOF
%       qNON  = nonlinear driving latent DOF
%       qNONrel = qNON / qLIN     (relative nonlinear amplitude)
%   The routine evaluates the precompiled regressor at qNONrel to obtain:
%       • w(qNONrel)    : current ECM weights
%       • d w / d q     : gradient of the weights wrt qNONrel
%   and then maps this 1D gradient to the pair (qLIN, qNON) via the chain rule.
%
% LATENT RATIO & SAFEGUARD
%   - If |qLIN| ≤ 1e−14, define qNONrel = 0 and set both partial derivatives
%     to zero to avoid division by ~0 and numerical blow-up.
%
% CHAIN-RULE MAPPING (for qLIN ≠ 0)
%   Let ẇ(q) = d w / d q at q = qNONrel. Then:
%     ∂w/∂qLIN = (∂qNONrel/∂qLIN) ẇ(qNONrel) = −(qNON/qLIN²) ẇ(qNONrel) = −(qNONrel/qLIN) ẇ(qNONrel)
%     ∂w/∂qNON = (∂qNONrel/∂qNON) ẇ(qNONrel) =  (1/qLIN)         ẇ(qNONrel)
%   These two vectors are returned (column-wise) in OPERFE.dw_q = [dw_lin, dw_non].
%
% INPUTS
%   VAR    : State container; uses VAR.DISP for the latent coordinates.
%   OPERFE : Must include wSTs_cluster.DATA_regress with fields:
%            • IndexLinear_damageMODEL     : position of qLIN in VAR.DISP
%            • IndexNonlinear_damageMODEL  : position of qNON in VAR.DISP
%            • nameFunctionEvaluate        : @(q, DATA_regress) → [w, dw_dq]
%
% OUTPUTS
%   OPERFE : Updated with
%            • wSTs  : ECM weights evaluated at qNONrel
%            • dw_q  : [∂w/∂qLIN, ∂w/∂qNON] (each column same size as w)
%   qNONrel: Scalar latent ratio used for the evaluation (0 if |qLIN| small)
%
% NOTES
%   • Dimensions: length(w) equals #ECM points; dw_lin and dw_non match that length.
%   • This helper assumes a **1-parameter regressor** in qNONrel; extend if the
%     weight law depends on additional latents.
%   • Downstream Jacobian builders (e.g., GetJacobianMatrix_MHROM) can include
%     the weight-sensitivity term only if enabled in the solver path.
%
% DEPENDENCIES
%   • OPERFE.wSTs_cluster.DATA_regress.nameFunctionEvaluate
%
% AUTHOR
%   Joaquín A. Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%
% COMMENTS UPDATED BY: ChatGPT-5
%   Date of update: 07-NOV-2025, Barcelona, Spain
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
end 

qLIN = VAR.DISP(OPERFE.wSTs_cluster.DATA_regress.IndexLinear_damageMODEL) ;   % first index
qNON =   VAR.DISP(OPERFE.wSTs_cluster.DATA_regress.IndexNonlinear_damageMODEL) ;  % Second index 

if abs(qLIN) <= 1e-14
    qNONrel = 0 ; 
else
    qNONrel = (qNON/qLIN) ; 
end

 
[w,dw_q]   = feval(OPERFE.wSTs_cluster.DATA_regress.nameFunctionEvaluate,qNONrel,OPERFE.wSTs_cluster.DATA_regress) ;
OPERFE.wSTs = w;


% \begin{equation}
%  \derpar{w_g(\qNONrel)}{\qLIN} = \derpar{\qNONrel}{\qLIN}   \dot{w}_g(\qNONrel) = - \dfrac{\qNONrel}{\qLIN} \dot{w}_g(\qNONrel)
% \end{equation}
% 
% \begin{equation}
%  \derpar{w_g(\qNONrel)}{\qNON} = \derpar{\qNONrel}{\qNON}   \dot{w}_g(\qNONrel) =   \dfrac{1}{\qLIN} \dot{w}_g(\qNONrel)
% \end{equation}

% if $\qLIN \neq 0$, otherwise both components are zero.

if abs(qLIN) <= 1e-14
    dw_lin = zeros(size(dw_q)) ; 
    dw_non = dw_lin ; 
else
     dw_lin =(-qNONrel/qLIN)*dw_q ; 
      dw_non= dw_q/qLIN ; 
end
 

OPERFE.dw_q = [dw_lin,dw_non];