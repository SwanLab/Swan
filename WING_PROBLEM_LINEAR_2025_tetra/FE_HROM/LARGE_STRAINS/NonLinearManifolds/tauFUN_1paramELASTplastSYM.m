function [tau,tauDER1,tauDER2] = tauFUN_1paramELASTplastSYM(qL,DATA)
%--------------------------------------------------------------------------
% tauFUN_1paramELASTplastSYM is a modification of tauFUN_1paramELASTplast,
% described below
% The modification is related to the fact that now the decoder is
% %
% \item Suppose we train the M-HROM in either tension and compression, and let us denote by $\signTRAIN$ the sign of $\qPL$ in this training set (so $\signTRAIN$ is either +1 or -1).
%
% \item The expression of the new tension/compression decoder reads
%
% \begin{equation}
% \begin{split}
%  \dL   =  \Decoder(\q) =   \BasisU{}{}  \tauNON(\q) & =  \BasisUel{}{} \qEL + \BasisUplMST{}{}\qPL +  \BasisUplSLV{}{}   \SIGN{\qPL} \signTRAIN \gPLslv(\absval{\qPL})
%  \\& =  \BasisUmst{}{} \q +  \BasisUplSLV{}{} \SIGN{\qPL} \signTRAIN \gPLslv(\absval{\qPL})
%  \end{split}
% \end{equation}
%
% Thus
%
% \begin{equation}
%  \tauNON(\q) = \coltres{\qEL}{\qPL}{\SIGN{\qPL} \signTRAIN \gPLslv(\absval{\qPL})}
% \end{equation}
% JAHO, 27-Oct-2025, Monday, Balmes 185, Barcelona
%.---------------------------------------------------


%
% PURPOSE of tauFUN_1paramELASTplast
%   Evaluate the ROM modal coefficients τ and their first/second derivatives
%   for a 2‑coordinate (elastic–plastic) manifold HROM. The reduced decoder
%   is split into:
%       • Elastic master coordinate          q_ELAST
%       • Plastic master coordinate          q_PLAST
%   The plastic slave coordinates are a learned nonlinear function of q_PLAST
%   represented via cubic B‑splines in a compact SVD basis.
%   Context: manifold-based ROMs and their hyperreduction (Secs. 12–13).

% MODEL (decoder in modal space)
%   Let qL = [ q_ELAST ; q_PLAST ].
%   The vector of decoder coefficients is
%       τ(qL) = [ q_ELAST ;
%                  q_PLAST ;
%                  U * ( S .* g(q_PLAST) ) ]
%   where g(·) collects the amplitudes of plastic slave modes in a compact
%   right-singular basis; U are left singular vectors and S the singular values
%   from the SVD of the slave snapshot matrix.
%
%   Derivatives (w.r.t. latent coordinates in the order [q_ELAST, q_PLAST])
%       τ₍qELAST₎  = [1 ; 0 ; 0]
%       τ₍qPLAST₎  = [0 ; 1 ; U * ( S .* g'(q_PLAST) )]
%       τ₍qELAST,qELAST₎ = 0,  τ₍qELAST,qPLAST₎ = 0,
%       τ₍qPLAST,qPLAST₎ = [0 ; 0 ; U * ( S .* g''(q_PLAST) )]
%
% INPUTS
%   qL    : [2×1] or [2×N] latent coordinates:
%             qL(1,:) = q_ELAST
%             qL(2,:) = q_PLAST
%   DATA  : struct produced by BsplinesLeastSquares_PLAST, containing:
%             .sp, .sp1, .sp2   : B‑form splines for g, g', g''
%             .xmin, .xmax      : spline training bounds for q_PLAST
%             .INFO_EXTRAP      : boundary values/derivs at xmin/xmax
%                                  (fields .XMIN.s0/s1/s2, .XMAX.s0/s1/s2)
%             .UleftSingular    : matrix U (slave left singular vectors)
%             .SSingular        : vector S (slave singular values)
%
% OUTPUTS
%   tau     : [n_modes × N] modal coefficients τ(qL)
%   tauDER1 : [n_modes × 2] (if N==1) or empty (if N>1 in this impl.)
%             First derivatives w.r.t. [q_ELAST, q_PLAST]
%   tauDER2 : [n_modes × 2 × 2] (if N==1) or empty (if N>1)
%             Second derivatives (Hessian entries per mode)
%
% METHOD
%   1) Evaluate g(q_PLAST), g'(q_PLAST), g''(q_PLAST) using
%      evaluate_spline_with_extrapolationFUNv with:
%         • cubic B‑splines in-domain,
%         • quadratic extrapolation outside [xmin, xmax] using DATA.INFO_EXTRAP.
%   2) Assemble τ and its derivatives using the formulas above, multiplying
%      the compact coefficients by U and S (Hadamard product with S).
%
% SHAPES & CONVENTIONS
%   • If qL is a column ([2×1]): returns τ, τ' and τ'' as full tensors
%     (sizes as stated above).
%   • If qL has multiple columns ([2×N]): returns τ for all N samples.
%     (For speed, τ' and τ'' are left empty in this vectorized branch.)
%
% ASSUMPTIONS / NOTES
%   • Small-strain elastoplastic setting; elastic contribution is linear in
%     q_ELAST; plastic slaves depend only on q_PLAST.
%   • DATA was previously trained with sorted, duplicate-free q_PLAST.
%   • Extrapolation is controlled by DATA.INFO_EXTRAP; verify that your
%     downstream solver tolerates queries beyond [xmin, xmax].
%
% PERFORMANCE
%   • Vectorized spline evaluation for all slave components; no per-mode loops.
%   • Compact SVD basis (U,S) minimizes the number of fitted scalar splines.
%
% EXAMPLE
%   % Given DATA from BsplinesLeastSquares_PLAST:
%   qL = [0.05 ; 0.10];
%   [tau, d1, d2] = tauFUN_1paramELASTplast(qL, DATA);
%
% VERSION
%   • 2025-08-13  (J.A. Hernández). Elastoplastic 2‑coordinate variant of
%     tauFUN_1paramVECT for manifold HROMs.
%--------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
    %     qL  = [-1;-0.3] ;
    %     qL  = [qL,qL] ;
end

indPLAST = DATA.IndexLatentPlasticVariable ;
signTRAIN = DATA.signTRAIN ;
signTRAIN = 1; 
if size(qL,2) == 1
    % Here qL = [qELAST; qPLAST_master]
    % We begin by evaluating qPLAST_slave as a function of qPLAST_master
    qPL_abs = abs(qL(indPLAST)) ;
   % qPL_abs = qL(indPLAST) ; % borrar esto 
%     if  qPL_abs < 0 
%         disp('Borra esto ')
%     end
    sign_qPL = sign(qL(indPLAST)) ;
    if sign_qPL == 0
        sign_qPL = 1 ;
    end
   % sign_qPL = 1;  % borrar esto 
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qPL_abs) ;
    
    nLATENT = size(qL,1) ;  % = 2
    %
    % \tauNON(\q) = \coltres{\qEL}{\qPL}{\SIGN{\qPL} \signTRAIN \gPLslv(\absval{\qPL})}
    
    gPLslv = DATA.UleftSingular*(DATA.SSingular.*f) ;
    gPLslvDER = DATA.UleftSingular*(DATA.SSingular.*df) ;
    gPLslvDERder = DATA.UleftSingular*(DATA.SSingular.*d2f) ;
    
    tau = [qL;  sign_qPL*signTRAIN*gPLslv] ;
    
    % First derivative
    % \begin{equation}
    %  \derpar{\tauNON(\q)}{\q} = \begin{bmatrix}
    %                                         1   &  0 \\
    %                                         0   &  1 \\
    %                                         0   &  \signTRAIN \gPLslvDER
    %                             \end{bmatrix}
    % \end{equation}
    
    
    tauDER1 = zeros(size(tau,1),nLATENT) ; % For instance, if there are 3 elastic modes, 1 nonlinear master modes  and 13 slave modes
    % the size is  17x4. Each row
    
    if nLATENT == 2
        % Version before 13-Oct-2025, only 1 elastic master variable
        iLATENT = 1;
        tauDER1(1,iLATENT) = 1;
        iLATENT = 2;
        tauDER1(2,iLATENT) = 1;
        % tauDER1(3:end,iLATENT) = signTRAIN*gPLslvDER;
        
        tauDER1(3:end,iLATENT) = signTRAIN*gPLslvDER;
        
        
        
        % NOW THE SECOND DERIVATIVE
        tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
        
        %
        %  \begin{equation}
        %  \begin{split}
        %   \tauNONderDER_i &      =  \SIGN{\qPL} \signTRAIN  \gPLslvDERder(\absval{\qPL})   \matcdos{0}
        %   {0}
        %   {0}
        %   {1}
        % \hspace{0.25cm}  i = 3,4 \ldots
        %  \end{split}
        %  \end{equation}
        
        
        
        tauDER2(3:end,2,2) = sign_qPL*signTRAIN*gPLslvDERder ;
        
        % tauDER2 = [zeros(size(qL));DATA.UleftSingular*(DATA.SSingular.*d2f)] ;
        
    else
        % FIRST DERIVATIVE S
        % General version
        % /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/14_HOM_PLAST.mlx
        % Recall that tau  are the amplitude of the modes in the linear
        % expansion
        % The first nLATENT entires are the amplitudes of the master modes
        % The remaining modes are the amplitudes of the slave modes
        % In this problem, all such slave coefficients are a function of
        % the indPLAST = nLATENT master mode   (the last master mode)
        % Accordingly, tauDER1 will adopt the expression
        %   tauDER1 = [der(tau)/der_q1, der(tau)/der_q2, ...der(tau)/der_q(indPLAST) ]
        % Thus
        tauDER1(1:nLATENT,1:nLATENT) = eye(nLATENT) ;
        tauDER1(nLATENT+1:end,indPLAST) =signTRAIN*gPLslvDER;
        
        
        
        % NOW THE SECOND DERIVATIVE
        %
        tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
        
        
        tauDER2(nLATENT+1:end,nLATENT,nLATENT) =sign_qPL*signTRAIN*gPLslvDERder ;;
        
        
        
    end
    
    
    
else
    % [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL(indPLAST,:)) ;
    
    
    qPL_abs = abs(qL(indPLAST,:)) ;
    sign_qPL = sign(qL(indPLAST,:));
    
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qPL_abs) ;
    
    gPLslv = DATA.UleftSingular*(DATA.SSingular.*f) ;
    
    sign_qPL_gPLslv = bsxfun(@times,gPLslv',sign_qPL' )' ;
    
    sign_qPL_gPLslv = gPLslv ; 
    tau = [qL;   signTRAIN*sign_qPL_gPLslv  ] ;
    
    tauDER1 = [] ;
    tauDER2 = [] ;
    
end



