function [tau,tauDER1,tauDER2] = tauFUN_1paramVECT_max(qL,DATA)
if nargin == 0
    load('tmp1.mat')
end
% tauFUN_1paramVECT_max is a modification of tauFUN_1paramVECT 
% The goal of the modification is to adapt it to the case in which there
% are no "master" modes
% % JAHO, 16-Nov-2025, Sunday, 10:52. Balmes 185, BArcelona
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/19_SHIFT_COMPRESS.mlx

%--------------------------------------------------------------------------
% tauFUN_1paramVECT
%
% PURPOSE
%   Computes the modal coefficients (tau) and their first and second
%   derivatives with respect to the single generalized coordinate qL
%   for a reduced-order model (ROM) with a nonlinear manifold.
%
%   The displacement field is expressed as:
%       d_L(qL) = Phi * tau(qL)
%               = [Phi_MASTER, Phi_SLAVE] * [ qL ; f(qL) ]
%   where:
%       - qL is the master coordinate (q_MASTER), scalar or vector.
%       - f(qL) is a nonlinear mapping for the slave coordinates,
%         parameterized via a spline fit.
%
%   The nonlinear function f(qL) is given by:
%       f(qL) = U * diag(S) * g(qL)
%   where g(qL) is evaluated via cubic B-form splines.
%
% INPUTS
%   qL                 : row/column vector of generalized coordinate values
%   DATA               : struct with fields:
%       .sp            : spline struct (B-form) for f(qL)
%       .sp1           : spline struct (B-form) for f'(qL)
%       .sp2           : spline struct (B-form) for f''(qL)
%       .xmin          : minimum valid spline abscissa
%       .xmax          : maximum valid spline abscissa
%       .INFO_EXTRAP   : struct with quadratic extrapolation data:
%                        .XMIN.s0, .XMIN.s1, .XMIN.s2, .XMAX.s0, .XMAX.s1, .XMAX.s2
%       .UleftSingular : matrix U from SVD (left singular vectors of slave basis)
%       .SSingular     : vector S from SVD (singular values for slave modes)
%
% OUTPUTS
%   tau      : modal coefficients tau(qL), size = [n_modes × length(qL)]
%   tauDER1  : first derivatives d(tau)/dqL, same size as tau
%   tauDER2  : second derivatives d²(tau)/dqL², same size as tau
%
% METHOD
%   1) Evaluate f(qL), f'(qL), f''(qL) using
%      evaluate_spline_with_extrapolationFUNv, which handles:
%       - In-domain cubic spline evaluation
%       - Quadratic extrapolation for qL outside [xmin, xmax]
%
%   2) Assemble tau and derivatives:
%       tau      = [ qL ;   U * (S .* f) ]
%       tauDER1  = [ 1  ;   U * (S .* f') ]
%       tauDER2  = [ 0  ;   U * (S .* f'') ]
%
% PERFORMANCE NOTES
%   - Uses the vectorized spline evaluation for all slave coordinates
%     in one call (f, df, d2f returned as matrices).
%   - Avoids looping over slave modes.
%
% EXAMPLE
%   [tau, dtau, d2tau] = tauFUN_1paramVECT(qL, DATA);
%
% AUTHOR / HISTORY
%   Joaquín A. Hernández, 10-Aug-2025, Molinos Marfagones (Cartagena)
%--------------------------------------------------------------------------
 
[f,df,d2f ] = evaluate_spline_with_extrapolationFUNv(DATA.sp, DATA.sp1, DATA.sp2, DATA.xmin,DATA.xmax,DATA.INFO_EXTRAP,qL) ;

 
tau = [DATA.UleftSingular*(DATA.SSingular.*f)] ;
tauDER1 = [DATA.UleftSingular*(DATA.SSingular.*df)] ;
tauDER2 = [DATA.UleftSingular*(DATA.SSingular.*d2f)] ;



