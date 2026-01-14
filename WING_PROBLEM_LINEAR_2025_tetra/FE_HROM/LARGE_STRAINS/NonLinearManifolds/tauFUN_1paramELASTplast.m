function [tau,tauDER1,tauDER2] = tauFUN_1paramELASTplast(qL,DATA)
if nargin == 0
    load('tmp1.mat')
    qL = qL(:,1) ;
    
end
% =========================================================================
% TAUFUN_1PARAMELASTPLAST — Decoder τ(q) and derivatives for EP manifold
% =========================================================================
% PURPOSE
%   Evaluate ROM modal coefficients τ and their first/second derivatives for a
%   2-coordinate (elastic–plastic) manifold HROM. The latent vector is
%       qL = [ q_ELAST ; q_PLAST ] ,
%   with plastic SLAVE coordinates represented as a nonlinear function of
%   q_PLAST via cubic B-splines fitted in a compact SVD basis.
%
% MODEL (decoder in modal space)
%   τ(qL) = [ q_ELAST ;
%             q_PLAST ;
%             U * ( S .* g(q_PLAST) ) ]
%   where g(·) are spline-fitted amplitudes of the plastic SLAVE coordinates
%   in the right-singular basis; U (left singular vectors) and S (singular
%   values) come from the compact SVD of SLAVE data.
%
%   First derivatives (w.r.t. [q_ELAST, q_PLAST]):
%     τ_{qELAST}  = [1 ; 0 ; 0]
%     τ_{qPLAST}  = [0 ; 1 ; U * ( S .* g'(q_PLAST) )]
%
%   Second derivatives:
%     τ_{qELAST,qELAST} = 0 ,  τ_{qELAST,qPLAST} = 0 ,
%     τ_{qPLAST,qPLAST} = [0 ; 0 ; U * ( S .* g''(q_PLAST) )]
%
% INPUTS
%   qL   : [nRED × 1] or [nRED × N] latent coordinates; for the standard
%          EP case nRED = 2 and rows are [q_ELAST; q_PLAST].
%   DATA : struct produced by BsplinesLeastSquares_PLAST (or SUBSAMPL) with:
%            .sp, .sp1, .sp2    — B-form splines for g, g′, g″
%            .xmin, .xmax       — spline training bounds in q_PLAST
%            .INFO_EXTRAP       — endpoint data for extrapolation
%                                  (.XMIN/.XMAX with fields s0/s1/s2)
%            .UleftSingular     — U (left singular vectors)
%            .SSingular         — S (singular values)
%            .IndexLatentPlasticVariable — index of q_PLAST within qL
%
% OUTPUTS
%   tau     : [n_modes × N] modal coefficients τ(qL).
%   tauDER1 : If N==1, [n_modes × nRED] first derivatives; empty otherwise.
%   tauDER2 : If N==1, [n_modes × nRED × nRED] second derivatives; empty otherwise.
%
% METHOD
%   1) Evaluate g, g′, g″ at q_PLAST using evaluate_spline_with_extrapolationFUNv:
%        • cubic B-splines inside [xmin, xmax],
%        • controlled quadratic extrapolation outside using DATA.INFO_EXTRAP.
%   2) Assemble τ, τ′, τ″ by multiplying compact coefficients with U and S
%      (Hadamard product with S).
%
% SHAPES & CONVENTIONS
%   • Column input (N=1): returns τ, τ′, τ″ fully populated.
%   • Multi-column input (N>1): returns τ for all samples; τ′, τ″ omitted
%     (vectorized fast path).
%
% ASSUMPTIONS / NOTES
%   • Small-strain elastoplastic setting; q_ELAST enters linearly; plastic
%     SLAVE coordinates depend only on q_PLAST.
%   • DATA must come from training with sorted/deduplicated q_PLAST.
%   • Extrapolation behavior should be verified for your online solver range.
%
% PERFORMANCE
%   • Vectorized evaluation for τ over N samples; compact SVD minimizes
%     the number of spline fits and multiplications.
%
% EXAMPLE
%   % DATA produced by BsplinesLeastSquares_PLAST:
%   qL = [0.05 ; 0.10];
%   [tau, d1, d2] = tauFUN_1paramELASTplast(qL, DATA);
%
% DEPENDENCIES
%   evaluate_spline_with_extrapolationFUNv (spline eval + extrapolation)
%
% VERSION HISTORY / AUTHORSHIP
%   • 2025-08-13 — Elastoplastic 2-coordinate variant of tauFUN_1paramVECT
%                  for manifold HROMs. (J.A. Hernández) — Molinos Marfagones, Cartagena.
%   • 2025-11-07 — Comments refreshed; clarified shapes, derivatives, and
%                  extrapolation details. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

%--------------------------------------------------------------------------


indPLAST = DATA.IndexLatentPlasticVariable ;

if size(qL,2) == 1
    % Here qL = [qELAST; qPLAST_master]
    % We begin by evaluating qPLAST_slave as a function of qPLAST_master
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL(indPLAST)) ;
    
    nLATENT = size(qL,1) ;  % = 2
    
    tau = [qL;DATA.UleftSingular*(DATA.SSingular.*f)] ;
    % Let us tackle now the first derivative. This is of the form
    %  tauDER1 = [d{tauDER1}/d{qELAST} , d{tauDER1}/d{qPLAST_master}     ],
    %  where
    % d{tauDER1}/d{qELAST} = [1,0,0....]^T
    % d{tauDER1}/d{qPLAST_master} = [0,1,UU*SS*df]^T
    tauDER1 = zeros(size(tau,1),nLATENT) ; % For instance, if there are 3 elastic modes, 1 nonlinear master modes  and 13 slave modes
    % the size is  17x4. Each row
    
    if nLATENT == 2
        % Version before 13-Oct-2025, only 1 elastic master variable
        iLATENT = 1;
        tauDER1(1,iLATENT) = 1;
        iLATENT = 2;
        tauDER1(2,iLATENT) = 1;
        tauDER1(3:end,iLATENT) = DATA.UleftSingular*(DATA.SSingular.*df);
        
        
        
        
        % NOW THE SECOND DERIVATIVE
        tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
        % where
        %tauDER2(1,:,:) = 0 ;
        %tauDER2(2,:,:) = 0 ;
        % tauDER2(j,:,:) = [d^2 tauDER_j/d{qELAST qELAST}   d^2 tauDER_j/d{qELAST qPLAST}
        %                   d^2 tauDER_j/d{qPLAST qELAST }  d^2 tauDER_j/d{qPLAST qPLAST}       ]
        %   =    [0  0
        %        0  d^2 tauDER_j/d{qPLAST qPLAST}       ]
        tauDER2(3:end,2,2) = DATA.UleftSingular*(DATA.SSingular.*d2f) ;
        
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
        tauDER1(nLATENT+1:end,indPLAST) = DATA.UleftSingular*(DATA.SSingular.*df);
        
        
        
        % NOW THE SECOND DERIVATIVE
        % 
        tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
        
        % where
        %tauDER2(1,:,:) = 0 ;
        %tauDER2(2,:,:) = 0 ;
        %tauDER2(nLATENT,:,:) = 0 ;
        %
        %
        % tauDER2(j,1,:) = [d^2{tau_j}/d{q1 q1}  &  d^2{tau_j}/d{q1 q2}   ....   d^2{tau_j}/d{q1 q_{nLATENT-1}}   d^2{tau_j}/d{q1 q_{nLATENT}}  }
        % tauDER2(j,2,:) = [d^2{tau_j}/d{q2 q1}  &  d^2{tau_j}/d{q2 q2}   ....   d^2{tau_j}/d{q2 q_{nLATENT-1}}   d^2{tau_j}/d{q2 q_{nLATENT}}  }    
        % .
        %.
        % . 
        % tauDER2(j,nL,:) = [d^2{tau_j}/d{q_nL q1}  &  d^2{tau_j}/d{q_nL q2}   ....   d^2{tau_j}/d{q_nL q_{nL-1}}   d^2{tau_j}/d{qnL q_NL}  }  
        
        
        tauDER2(nLATENT+1:end,nLATENT,nLATENT) = DATA.UleftSingular*(DATA.SSingular.*d2f) ;
        
        
        
    end
    
    
    
else
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL(indPLAST,:)) ;
    
    
    tau = [qL;DATA.UleftSingular*(DATA.SSingular.*f)] ;
    
    tauDER1 = [] ;
    tauDER2 = [] ;
    
end



