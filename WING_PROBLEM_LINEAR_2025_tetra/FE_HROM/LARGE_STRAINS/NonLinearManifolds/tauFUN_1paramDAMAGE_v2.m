function [tau,tauDER1,tauDER2] = tauFUN_1paramDAMAGE_v2(qL,DATA)
% tauFUN_1paramDAMAGE_v2
% -------------------------------------------------------------------------
% PURPOSE
%   Evaluate the reduced coefficients τ(q) and their first/second derivatives
%   for a *two-latent-variable* (qLIN, qNON) manifold ROM in the univariate
%   damage setting. This version modifies the original `tauFUN_1paramDAMAGE`
%   to make the differentiation well-posed near qLIN ≈ 0 by reformulating
%   the decoder in terms of the *relative coordinate* qNONrel = qNON/qLIN
%   and introducing tolerance-based regularization.
%
% THEORETICAL BACKGROUND (document pp. 138–146)
%   • Decoder structure:
%         τ(q) = [ qLIN ;
%                   qNON ;
%                   qLIN * g(qNON/qLIN) ],
%     where g(·) stacks the slave-mode amplitudes represented in a compact
%     SVD basis (U,S). The physical displacements are reconstructed as
%         d ≈ Φ τ(q).
%
%   • Key idea of v2 formulation:
%       The original decoder τ(q) = [qLIN; qLIN*qNON; qLIN*g(qNON)] suffers
%       from undefined derivatives at qLIN=0. By expressing the nonlinear
%       term as g(qNON/qLIN), all derivatives remain finite and continuous
%       (except at qLIN=0, where they are set to 0 following the true
%       function behavior τ_i=0 when qLIN=0).
%
%   • Derivatives (columns ordered as [qLIN, qNON]):
%
%       First derivative:
%         ∂τ/∂q =
%           [ 1   0 ;
%             0   1 ;
%             g(qNONrel) - qNONrel*g'(qNONrel)    g'(qNONrel) ] ,
%
%       Second derivative (per slave component i≥3):
%         ∂²τ_i/∂q² =
%           [ (qNONrel²/qLIN) g''(qNONrel) ,  -(qNONrel/qLIN) g''(qNONrel) ;
%             -(qNONrel/qLIN) g''(qNONrel) ,   (1/qLIN) g''(qNONrel) ] .
%
%       For qLIN = 0, all second-derivative entries are set to 0,
%       consistent with τ_i(qLIN=0,qNON)=0.
%
% INPUTS
%   qL   : [2×1] or [2×N] vector/matrix of latent coordinates
%           qL(1,:) = qLIN  (linear/elastic coordinate)
%           qL(2,:) = qNON  (nonlinear/damage coordinate)
%
%   DATA : Struct containing spline and basis data:
%           .sp, .sp1, .sp2   → cubic B-splines for g, g′, g′′
%           .xmin, .xmax      → spline training bounds
%           .INFO_EXTRAP      → extrapolation data (values & derivs)
%           .UleftSingular    → left singular vectors (U)
%           .SSingular        → singular values (S)
%
% OUTPUTS
%   tau     : [nModes×N] decoder coefficients τ(q).
%   tauDER1 : [nModes×2] first derivatives wrt [qLIN, qNON] (only for N=1).
%   tauDER2 : [nModes×2×2] second derivatives per mode (only for N=1).
%
% METHOD
%   1) Compute the *relative coordinate* qNONrel = qNON/qLIN.
%      If |qLIN| < TOL, set qNONrel = 0 (regularization).
%   2) Evaluate the spline-based nonlinear function and its derivatives:
%        [f, df, d2f] = evaluate_spline_with_extrapolationFUNv(..., qNONrel).
%   3) Build the slave term: gNONslv = U * (S .* f).
%   4) Assemble τ(q):
%         τ = [qLIN ; qNON ; qLIN*gNONslv].
%   5) Compute first derivatives:
%         τ_qLIN  = [1; 0; g - qNONrel*g′],
%         τ_qNON  = [0; 1; g′].
%   6) Compute second derivatives (Hessian entries) for slave components:
%         τ_{11} = (qNONrel²/qLIN) g′′,
%         τ_{12} = τ_{21} = -(qNONrel/qLIN) g′′,
%         τ_{22} = (1/qLIN) g′′,
%      setting all to 0 if |qLIN| < TOL.
%
% NUMERICAL NOTES
%   • The use of qNONrel ensures smooth derivatives everywhere except at
%     qLIN=0, where the function is exactly 0 (no discontinuity).
%   • Setting Hessian terms to 0 at qLIN=0 matches the true τ(q) behavior.
%   • Compact SVD basis (U,S) avoids per-mode spline evaluations.
%
% VERSION / AUTHOR
%   24-Oct-2025 — J.A. Hernández (Barcelona).
%   Revised decoder τ(q) = [qLIN; qNON; qLIN*g(qNON/qLIN)] with corrected
%   chain-rule implementation and safe qLIN→0 treatment.
% -------------------------------------------------------------------------
% -----------------------------------------
if nargin == 0
    load('tmp1.mat')
    % qL = qL(:,1) ;
    
    
end

%indNON = DATA.IndexLatentNonlinearVariable ;  % iT SHOULD BE =2
TOL_select= 1e-12 ;

if size(qL,2) == 1
    
    % tau here are the amplitudes of the modal basis
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/MLEARNstruct_1.pdf
    
    %  \tauNON(\q) \defeq  \coltres{\qLIN}{  \qNON}{\qLIN \gNONslv(\qNONrel)}
    % where
    %\qNONrel = \dfrac{\qNON}{\qLIN}
    %  $\qLIN \neq 0$  (otherwise $\qNONrel = 0$)
    
    iLINEAR= 1;
    qLIN = qL(iLINEAR,:) ;
    % Now the second component is
    indNON = 2;
    qNON = qL(indNON,:) ;
    % Now we have to compute qNONrel = qNON/qLIN. WE want to avoid
    % ill-posed cases, in which qLIN \approx 0
    if abs(qLIN)<= TOL_select
        qNONrel = 0 ;
    else
        qNONrel = qNON/qLIN ;
    end
    
    % We begin by evaluating qPLAST_slave as a function of qPLAST_master
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qNONrel) ;
    
    nLATENT = size(qL,1) ;
    % Finally
    gNONslv = DATA.UleftSingular*(DATA.SSingular.*f) ;
    % \tauNON(\q) \defeq  \coltres{\qLIN}{  \qNON}{\qLIN \gNONslv(\qNONrel)}
    tau = [qLIN;
        qNON ;
        qLIN*gNONslv] ;
    
    
    % FIRST DERIVATIVE
    % ***********************************************+
    %
    %     \begin{equation}
    % \label{eq:24Oct2025:cuba}
    % \begin{split}
    %   \tauNONder &  =
    %           \\& =   \begin{bmatrix}
    %                             1  &  0  \\
    %                             0  &    1  \\
    %                              \gNONslv(\qNONrel) - \qNONrel \gNONslvDER(\qNONrel)   &   \gNONslvDER(\qNONrel)
    %                   \end{bmatrix}
    %  \end{split}
    % \end{equation}
    
    %
    gNONslvDER = DATA.UleftSingular*(DATA.SSingular.*df);
    
    tauDER1 = zeros(size(tau,1),nLATENT) ;
    %  1  &  0  \\
    tauDER1(1,:) = [1,0] ;
    % 0  &    1  \\
    tauDER1(2,:)  = [0,1] ;
    %   \gNONslv(\qNONrel) - \qNONrel \gNONslvDER(\qNONrel)   &   \gNONslvDER(\qNONrel)
    tauDER1(3:end,:) = [gNONslv - qNONrel*gNONslvDER, gNONslvDER] ;
    
    
    
    % SECOND DERIVATIVE
    % ...................
    gNONslvDERder = DATA.UleftSingular*(DATA.SSingular.*d2f);
    tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
    
    %  iCOMP = 1; % First component of tau ...ALL zero,
    % iCOMP = 2; % Second component of tau ...ALL zero
    
    % Finally, remaining entries
    iCOMP = 3:size(tau,1);
    % *********************************************************************
    % Component 1,1 =  \dfrac{\qNONrel^2}{\qLIN} \gNONslvDERder(\qNONrel)
    % *********************************************************************+
    if abs(qLIN)<= TOL_select
        % NOTHING IS DONE, SET TO ZERO
    else
        tauDER2(iCOMP,1,1) = (qNONrel^2/qLIN)*gNONslvDERder ;
    end
    % *********************************************************************
    % Component 1,2 and 2,1 =  - \dfrac{\qNONrel}{\qLIN}       \gNONslvDERder(\qNONrel)
    % *********************************************************************+
    if abs(qLIN)<= TOL_select
        % NOTHING IS DONE, SET TO ZERO
    else
        tauDER2(iCOMP,1,2) =  (-qNONrel/qLIN)*gNONslvDERder ;
        tauDER2(iCOMP,2,1) =  tauDER2(iCOMP,1,2) ;
    end
    % *********************************************************************
    % Component 2,2 =  \dfrac{\gNONslvDERder(\qNONrel)}{\qLIN}
    % *********************************************************************+
    if abs(qLIN)<= TOL_select
        % NOTHING IS DONE, SET TO ZERO
    else
        tauDER2(iCOMP,2,2) =  gNONslvDERder/qLIN ;
    end
    
    
    
    
else
    
    
    iLINEAR= 1;
    qLIN = qL(iLINEAR,:) ;
    % Now the second component is
    indNON = 2;
    qNON = qL(indNON,:) ;
    
    qNONrel  = zeros(size(qNON)) ;
    IndicesNonZero = find(abs(qLIN) >= TOL_select )  ;
    qNONrel(IndicesNonZero) =  qNON(IndicesNonZero)./qLIN(IndicesNonZero) ;
    
    PLOT = 0; 
    if PLOT == 1
        figure(4566)
        hold on
        xlabel('Step')
        ylabel('q')
        plot(qLIN,'DisplayName','qLIN')
        plot(qNON,'DisplayName','qNON')
        plot(qNONrel,'DisplayName','qNONrel')
        legend show
    end
    
    
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qNONrel) ;
    
    
    gNONslv = DATA.UleftSingular*(DATA.SSingular.*f) ;
    
    %  \tauNON(\q) \defeq  \coltres{\qLIN}{  \qNON}{\qLIN \gNONslv(\qNONrel)}
    % where
    %\qNONrel = \dfrac{\qNON}{\qLIN}
    %  $\qLIN \neq 0$  (otherwise $\qNONrel = 0$)
    
    tau = zeros(2+size(gNONslv,1),size(qL,2)) ;
    tau(1,:) = qLIN  ;
    tau(2,:) =  qNON  ;
    qLINrep = repmat(qLIN,size(gNONslv,1),1);
    tau(3:end,:) =  qLINrep.*gNONslv ;
    
    
    %
    %     tau = [qL;DATA.UleftSingular*(DATA.SSingular.*f)] ;
    
    tauDER1 = [] ;
    tauDER2 = [] ;
    
end



