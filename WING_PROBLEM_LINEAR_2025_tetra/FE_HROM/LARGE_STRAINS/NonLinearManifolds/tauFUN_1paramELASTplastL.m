function [tau,tauDER1,tauDER2] = tauFUN_1paramELASTplast(qL,DATA)
if nargin == 0
    load('tmp1.mat')
    qL = [-0.2;0.2   ] ;
end
%--------------------------------------------------------------------------
% tauFUN_1paramELASTplast
%
% PURPOSE
%   Evaluate the reduced modal coefficients τ and their first/second
%   derivatives for a 2-coordinate manifold ROM where the decoder (in
%   the displacement space) is:
%
%       d_L(q_ELAST, q_PLAST) = PhiElast * q_ELAST + PhiNON * g(q_PLAST)
%
%   This function returns τ(q) in the modal space that pairs with the
%   decoder above:
%
%       τ(q) = [ q_ELAST ;
%                U * ( S .* g(q_PLAST) ) ]
%
%   i.e., the elastic contribution is linear in a single master coordinate
%   q_ELAST, while the plastic contribution is a learned nonlinear map of
%   q_PLAST expressed in a compact SVD basis (U, S).
%
% INPUTS
%   qL    : [2×1] or [2×N] array of latent coordinates
%           qL(1,:) = q_ELAST
%           qL(2,:) = q_PLAST
%
%   DATA  : struct with fields produced during spline training for the
%           plastic map g(·):
%             • sp, sp1, sp2   : B-form splines for g, g', g''
%             • xmin, xmax     : training bounds for q_PLAST
%             • INFO_EXTRAP    : boundary data for controlled extrapolation
%                                 (fields XMIN.s0/s1/s2 and XMAX.s0/s1/s2)
%             • UleftSingular  : matrix U (left singular vectors of slaves)
%             • SSingular      : vector S (singular values of slaves)
%
% OUTPUTS
%   tau     : [n_modes × N] modal coefficients τ(q)
%
%   tauDER1 : First derivatives w.r.t. [q_ELAST, q_PLAST]
%             • If N == 1 → [n_modes × 2] with
%                 dτ/dq_ELAST = [1 ; 0_(n_modes-1)]
%                 dτ/dq_PLAST = [0 ; U * (S .* g'(q_PLAST))]
%             • If N  > 1 → [] (omitted for speed in the vectorized path)
%
%   tauDER2 : Second derivatives (only along q_PLAST in this model)
%             • If N == 1 → [n_modes × 2 × 2] with
%                 ∂²τ/∂q_ELAST²   = 0
%                 ∂²τ/∂q_ELAST∂q_PLAST = 0
%                 ∂²τ/∂q_PLAST²  = [0 ; U * (S .* g''(q_PLAST))]
%             • If N  > 1 → [] (omitted for speed)
%
% METHOD
%   1) Evaluate g, g', g'' at q_PLAST using
%        evaluate_spline_with_extrapolationFUNv(sp, sp1, sp2, xmin, xmax, INFO_EXTRAP, ·)
%      which applies cubic B-splines in-domain and controlled extrapolation
%      outside [xmin, xmax] using INFO_EXTRAP.
%   2) Form τ, dτ/dq, and d²τ/dq² by multiplying the compact coefficients
%      with U and applying element-wise scaling by S.
%
% SHAPE NOTES
%   • n_modes = 1 (elastic master) + size(U,1) (plastic slave block).
%   • When qL has multiple columns (N>1), τ is returned for all samples,
%     while derivatives are skipped for efficiency.
%
% ASSUMPTIONS / SCOPE
%   • Small-strain elastoplastic setting; only the plastic block is nonlinear.
%   • The decoder does not include a linear q_PLAST term (slaves only).
%   • Extrapolation behavior is governed by DATA.INFO_EXTRAP; ensure your
%     solver tolerates calls beyond the training range.
%
% PERFORMANCE
%   • Fully vectorized spline evaluation over slave components.
%   • Compact SVD basis (U,S) minimizes the number of fitted splines.
%
% VERSION
%   • 2025-08-18 — Updated to decoder d_L = PhiElast*q_ELAST + PhiNON*g(q_PLAST)
%                  (no explicit linear q_PLAST term in τ).
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



if size(qL,2) == 1
    % Here qL = [qELAST; qPLAST_master]
    % We begin by evaluating qPLAST_slave as a function of qPLAST_master
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL(2)) ;
    
    nLATENT = size(qL,1) ;  % = 2
    
    tau = [qL(1);DATA.UleftSingular*(DATA.SSingular.*f)] ;
    % Let us tackle now the first derivative. This is of the form
    %  tauDER1 = [d{tauDER1}/d{qELAST} , d{tauDER1}/d{qPLAST_master}     ],
    %  where
    % d{tauDER1}/d{qELAST} = [1,0,0....]^T
    % d{tauDER1}/d{qPLAST_master} = [0,1,UU*SS*df]^T
    tauDER1 = zeros(size(tau,1),nLATENT) ;
    iLATENT = 1;
    tauDER1(1,iLATENT) = 1;
    iLATENT = 2;
 %   tauDER1(2,iLATENT) = 1;
    tauDER1(2:end,iLATENT) = DATA.UleftSingular*(DATA.SSingular.*df);
    
    % NOW THE SECOND DERIVATIVE
    tauDER2 = zeros(size(tau,1),nLATENT,nLATENT) ;
    % where
    %tauDER2(1,:,:) = 0 ;
    %tauDER2(2,:,:) = 0 ;
    % tauDER2(j,:,:) = [d^2 tauDER_j/d{qELAST qELAST}   d^2 tauDER_j/d{qELAST qPLAST}
    %                   d^2 tauDER_j/d{qPLAST qELAST }  d^2 tauDER_j/d{qPLAST qPLAST}       ]
    %   =    [0  0
    %        0  d^2 tauDER_j/d{qPLAST qPLAST}       ]
    tauDER2(2:end,2,2) = DATA.UleftSingular*(DATA.SSingular.*d2f) ;
    
    % tauDER2 = [zeros(size(qL));DATA.UleftSingular*(DATA.SSingular.*d2f)] ;
    
else
    [f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL(2,:)) ;
    
    
    tau = [qL(1,:);DATA.UleftSingular*(DATA.SSingular.*f)] ;
    
    tauDER1 = [] ;
    tauDER2 = [] ;
    
end



