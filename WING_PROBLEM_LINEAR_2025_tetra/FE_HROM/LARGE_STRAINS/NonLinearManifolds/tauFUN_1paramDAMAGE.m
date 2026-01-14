function [tau,tauDER1,tauDER2] = tauFUN_1paramDAMAGE(qL,DATA)
% tauFUN_1paramDAMAGE
%--------------------------------------------------------------------------
% PURPOSE
%   Evaluate the modal coefficients τ(q) and their first/second derivatives
%   for a 2-latent-coordinate manifold ROM tailored to *univariate damage*
%   problems at small strains. The decoder uses one “linear/elastic” master
%   coordinate (qLIN) and one “nonlinear/structural” master coordinate
%   (qNON); the remaining slave amplitudes are nonlinear functions of qNON
%   represented in a compact SVD basis.
%
% THEORETICAL BACKGROUND (see *MLEARNstruct_1.pdf*, pp. 138–146)
%   • Decoder structure for univariate damage:
%         τ(q) = [ qLIN ;
%                   qLIN*qNON ;
%                   qLIN * g(qNON) ] ,
%     where g(·) stacks the slave-mode amplitudes in the right-singular
%     basis; the physical displacements are d ≈ Φ τ(q). :contentReference[oaicite:0]{index=0}
%
%   • First derivatives (columns ordered as [qLIN, qNON]):
%         ∂τ/∂q = [ 1        0 ;
%                   qNON   qLIN ;
%                   g(qNON)  qLIN*g'(qNON) ] .  :contentReference[oaicite:1]{index=1}
%
%   • Second derivatives (Hessian blocks per component τi):
%       – For τ1 = qLIN:          [0 0; 0 0].                                  :contentReference[oaicite:2]{index=2}
%       – For τ2 = qLIN*qNON:     [0 1; 1 0].                                  :contentReference[oaicite:3]{index=3}
%       – For τi≥3 = qLIN*g_i(qNON):
%                                 [0 g'_i; g'_i  qLIN*g''_i].                  :contentReference[oaicite:4]{index=4}
%     These τ̈i blocks enter the geometric term of the tangent stiffness in
%     the manifold formulation. :contentReference[oaicite:5]{index=5}
%
% INPUTS
%   qL   : [2×1] or [2×N] latent coordinates with
%           qL(1,:) = qLIN ,  qL(2,:) = qNON .
%   DATA : struct with fields
%           .sp, .sp1, .sp2     — cubic B-spline fits of g, g', g''
%           .xmin, .xmax        — training bounds for qNON
%           .INFO_EXTRAP        — boundary values/derivs for safe extrapolation
%           .UleftSingular      — left singular vectors (U) for slave subspace
%           .SSingular          — singular values (S) matching U
%           .IndexLatentNonlinearVariable — index of qNON (expected = 2)
%
% OUTPUTS
%   tau     : [nModes × N]   modal coefficients τ(q).
%   tauDER1 : [nModes × 2]   first derivatives w.r.t. [qLIN, qNON] (N=1 case).
%   tauDER2 : [nModes × 2 × 2] second derivatives per mode (N=1 case).
%             (For vectorized N>1 evaluation, tauDER1/tauDER2 are omitted.)
%
% ASSUMPTIONS & NOTES
%   • Small-strain damage setting; two latent coordinates suffice to capture
%     the linear response and a single “structural” nonlinear evolution
%     coordinate, per the univariate damage ROM construction. :contentReference[oaicite:6]{index=6}
%   • Decoder/encoder orthogonality uses a metric G (often K_lin,ll) for
%     modal normalization and projection in the manifold ROM. :contentReference[oaicite:7]{index=7}
%   • Slave amplitudes g(qNON) are evaluated via spline fits with controlled
%     extrapolation outside [xmin, xmax]; derivatives follow from the fitted
%     splines (g', g'').
%
% PERFORMANCE
%   • Evaluations are fully vectorized over slave coordinates using the
%     compact (U,S) representation to avoid per-mode spline calls.
%
% VERSION / AUTHOR
%   • 24-Oct-2025 — Adapted from the elastoplastic τ-map to the univariate
%     damage case (J.A. Hernández). See cited pages for formulas. :contentReference[oaicite:8]{index=8}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
   % qL = qL(:,1) ;
   qL = [0,0]' ; 
    
end

indNON = DATA.IndexLatentNonlinearVariable ;  % iT SHOULD BE =2

if size(qL,2) == 1
    % Here qL = [qELAST; qPLAST_master]
    % We begin by evaluating qPLAST_slave as a function of qPLAST_master
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL(indNON)) ;
    
    nLATENT = size(qL,1) ;  % = 2
    
    % tau here are the amplitudes of the modal basis
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/DEVELOPMENTS/MLEARNstruct_1.pdf
    
    %
    %     \begin{equation}
    % \label{eq:23oct2025_npm}
    %  \tauNON(\q) \defeq  \coltres{\qLIN}{\qLIN \qNON}{\qLIN \gNONslv(\qNON)}
    % \end{equation}
    iLINEAR= 1;
    qLIN = qL(iLINEAR,:) ;
    % Now the second component is
    qNON = qL(indNON,:) ;
    % Finally
    gNONslv = DATA.UleftSingular*(DATA.SSingular.*f) ;
    tau = [qLIN;
        qLIN*qNON ;
        qLIN*gNONslv] ;
    
    %%%%% FIRST DERIVATIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    \begin{equation}
    %   \tauNONder  = \derpar{\tauNON}{\q} =    \begin{bmatrix}
    %                             \derpar{\qLIN}{\qLIN}  &  \derpar{\qLIN}{\qNON}  \\
    %                             \derpar{\qLIN \qNON}{\qLIN}  &  \derpar{\qLIN \qNON}{\qNON}  \\
    %                             \derpar{\qLIN \gNONslv(\qNON)}{\qLIN}  &  \derpar{\qLIN \gNONslv(\qNON)}{\qNON}  \\
    %                   \end{bmatrix}  =
    % \begin{bmatrix}
    %                             1  &  0  \\
    %                             \qNON  &  \qLIN  \\
    %                              \gNONslv(\qNON)   &  \qLIN  \gNONslvDER(\qNON)  \\
    %                   \end{bmatrix}_{\nmodesU \times \nRED}
    % \end{equation}
    
    %
    gNONslvDER = DATA.UleftSingular*(DATA.SSingular.*df);
    gNONslvDERder = DATA.UleftSingular*(DATA.SSingular.*d2f);
    
    tauDER1 = zeros(size(tau,1),nLATENT) ;
    %        1  &  0  \\
    % %     \qNON  &  \qLIN  \\
    % %    \gNONslv(\qNON)   &  \qLIN  \gNONslvDER(\qNON)  \\
    tauDER1(:,1) = [1;qNON; gNONslvDER];
    tauDER1(:,2) = [0;qLIN; qLIN*gNONslvDER];
    
    % SECOND DERIVATIVE
    tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
    
    %  iCOMP = 1; % First component of tau ...ALL zero, see below
    %
    %  \begin{equation}
    %   \tauNONderDER_1  = \matcdos{\dderparIND{\qLIN}{\qLIN}{\qLIN}}
    %  {\dderparIND{\qLIN}{\qLIN}{\qNON}}
    %   {\dderparIND{\qLIN}{\qLIN}{\qNON}}
    %  {\dderparIND{\qLIN}{\qNON}{\qNON}}  = \matcdos{0}{0}{0}{0}
    %  \end{equation}
    iCOMP = 2; % Second component of tau
    %
    %  \begin{equation}
    %   \tauNONderDER_2  = \matcdos{\dderparIND{\qLIN \qNON}{\qLIN}{\qLIN}}
    %  {\dderparIND{\qLIN\qNON}{\qLIN}{\qNON}}
    %   {\dderparIND{\qLIN\qNON}{\qLIN}{\qNON}}
    %  {\dderparIND{\qLIN\qNON}{\qNON}{\qNON}}  = \matcdos{0}{1}{1}{0}
    %  \end{equation}
    
    tauDER2(iCOMP,:,:) = [0  1; 1 0] ;
    
    % Finally, remaining entrie
    iCOMP = 3:size(tau,1) ;
    %     \begin{equation}
    %   \tauNONderDER_i  = \matcdos{\dderparIND{\qLIN \gNONslv_{i-2}(\qNON)}{\qLIN}{\qLIN}}
    %  {\dderparIND{\qLIN\gNONslv_{i-2}(\qNON)}{\qLIN}{\qNON}}
    %   {\dderparIND{\qLIN\gNONslv_{i-2}(\qNON)}{\qLIN}{\qNON}}
    %  {\dderparIND{\qLIN\gNONslv_{i-2}(\qNON)}{\qNON}{\qNON}}  =
    
    %\matcdos{0}
    %{\gNONslvDER_{i-2}(\qNON)}
    %{\gNONslvDER_{i-2}(\qNON)}
    % {\qLIN \gNONslvDERder_{i-2}(\qNON) }, \hspace{0.25cm}  i = 3,4 \ldots
    %  \end{equation}
    
    tauDER2(iCOMP,1,2) = gNONslvDER ;
    tauDER2(iCOMP,2,1) = gNONslvDER ;
    tauDER2(iCOMP,2,2) = qLIN*gNONslvDERder ;
    
    
    
    %     % where
    %     %tauDER2(1,:,:) = 0 ;
    %     %tauDER2(2,:,:) = 0 ;
    %     % tauDER2(j,:,:) = [d^2 tauDER_j/d{qELAST qELAST}   d^2 tauDER_j/d{qELAST qPLAST}
    %     %                   d^2 tauDER_j/d{qPLAST qELAST }  d^2 tauDER_j/d{qPLAST qPLAST}       ]
    %     %   =    [0  0
    %     %        0  d^2 tauDER_j/d{qPLAST qPLAST}       ]
    %     tauDER2(3:end,2,2) = DATA.UleftSingular*(DATA.SSingular.*d2f) ;
    
    
    %  if nLATENT == 2
    % Version before 13-Oct-2025, only 1 elastic master variable
    
    
    % tauDER2 = [zeros(size(qL));DATA.UleftSingular*(DATA.SSingular.*d2f)] ;
    
    %     else
    %         % FIRST DERIVATIVE S
    %         % General version
    %         % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/14_HOM_PLAST.mlx
    %         % Recall that tau  are the amplitude of the modes in the linear
    %         % expansion
    %         % The first nLATENT entires are the amplitudes of the master modes
    %         % The remaining modes are the amplitudes of the slave modes
    %         % In this problem, all such slave coefficients are a function of
    %         % the indPLAST = nLATENT master mode   (the last master mode)
    %         % Accordingly, tauDER1 will adopt the expression
    %         %   tauDER1 = [der(tau)/der_q1, der(tau)/der_q2, ...der(tau)/der_q(indPLAST) ]
    %         % Thus
    %         tauDER1(1:nLATENT,1:nLATENT) = eye(nLATENT) ;
    %         tauDER1(nLATENT+1:end,indPLAST) = DATA.UleftSingular*(DATA.SSingular.*df);
    %
    %
    %
    %         % NOW THE SECOND DERIVATIVE
    %         %
    %         tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
    %
    %         % where
    %         %tauDER2(1,:,:) = 0 ;
    %         %tauDER2(2,:,:) = 0 ;
    %         %tauDER2(nLATENT,:,:) = 0 ;
    %         %
    %         %
    %         % tauDER2(j,1,:) = [d^2{tau_j}/d{q1 q1}  &  d^2{tau_j}/d{q1 q2}   ....   d^2{tau_j}/d{q1 q_{nLATENT-1}}   d^2{tau_j}/d{q1 q_{nLATENT}}  }
    %         % tauDER2(j,2,:) = [d^2{tau_j}/d{q2 q1}  &  d^2{tau_j}/d{q2 q2}   ....   d^2{tau_j}/d{q2 q_{nLATENT-1}}   d^2{tau_j}/d{q2 q_{nLATENT}}  }
    %         % .
    %         %.
    %         % .
    %         % tauDER2(j,nL,:) = [d^2{tau_j}/d{q_nL q1}  &  d^2{tau_j}/d{q_nL q2}   ....   d^2{tau_j}/d{q_nL q_{nL-1}}   d^2{tau_j}/d{qnL q_NL}  }
    %
    %
    %         tauDER2(nLATENT+1:end,nLATENT,nLATENT) = DATA.UleftSingular*(DATA.SSingular.*d2f) ;
    %
    %
    %
    %     end
    %
    
    
else
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL(indNON,:)) ;
    
    
    iLINEAR= 1;
    qLIN = qL(iLINEAR,:) ;
    % Now the second component is
    qNON = qL(indNON,:) ;
    % Finally
    gNONslv = DATA.UleftSingular*(DATA.SSingular.*f) ;
    
    tau = zeros(2+size(gNONslv,1),size(qL,2)) ; 
    tau(1,:) = qLIN  ; 
    tau(2,:) = qLIN.*qNON  ; 
    qLINrep = repmat(qLIN,size(gNONslv,1),1);
    tau(3:end,:) =  qLINrep.*gNONslv ; 
    
 
    %
    %     tau = [qL;DATA.UleftSingular*(DATA.SSingular.*f)] ;
    
    tauDER1 = [] ;
    tauDER2 = [] ;
    
end



