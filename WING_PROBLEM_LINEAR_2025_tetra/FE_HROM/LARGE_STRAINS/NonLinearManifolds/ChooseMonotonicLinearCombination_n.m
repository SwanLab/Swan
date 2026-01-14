function   [y,c] = ChooseMonotonicLinearCombination_n(Q,n)
%CHOOSEMONOTONICLINEARCOMBINATION Smoothest/most-monotone combination in span(Q).
%   y = CHOOSEMONOTONICLINEARCOMBINATION(Q) returns the vector y = Q*c that
%   minimizes the discrete first-difference energy ||D*Q*c||_2^2 subject to
%   ||Q*c||_2^2 = 1, where D is the (M-1)-by-M forward-difference matrix.
%   This is a Rayleigh–Ritz / generalized eigenvalue problem:
%       minimize   c' * (Q'*D'*D*Q) * c   /   c' * (Q'*Q) * c
%     whose minimizer is the smallest generalized eigenvector of
%       (Q'*D'*D*Q) c = lambda (Q'*Q) c.
%   The resulting y typically exhibits the smoothest trend within range(Q)
%   and, in many practical cases, is (nearly) monotone.
%
% INPUT
%   Q : (M-by-N) matrix whose columns are sampled basis functions q_j(F_i)
%       stacked over M grid points F_1 < ... < F_M. Each column is treated
%       as a length-M vector; y will be a linear combination of these.
%
% OUTPUT
%   y : (M-by-1) vector y = Q*c, the smoothest (often most monotone) vector
%       in the span of Q under the discrete first-difference metric.
%
% THEORY (brief)
%   Let D be the forward-difference operator so that (D*y)_k = y_{k+1}-y_k.
%   The quadratic form y'*(D'*D)*y = ||D*y||_2^2 penalizes oscillations.
%   Restricting y to the subspace range(Q) and normalizing ||y||_2 = 1
%   yields a Rayleigh quotient whose minimizer is found via the generalized
%   eigenproblem (A,B) with A = Q'*D'*D*Q and B = Q'*Q.
%
% PRACTICAL NOTES
%   • Grid spacing: if the F-grid is nonuniform, replace D by Dh
%     (divided differences) with Dh = diag(1./diff(F)) * D to account for
%     variable step sizes.
%   • Numerical rank: if Q is ill-conditioned or rank-deficient, consider
%     orthonormalizing first (e.g., [Z,~] = qr(Q,0); then solve on Z).
%   • Direction: the solution is defined up to sign. If you prefer strictly
%     increasing behavior, you may flip y so that sum(D*y) >= 0.
%   • Exact monotonicity: this method *encourages* monotonicity. To *enforce*
%     Dy >= 0, solve a small convex QP/SOCP with linear inequalities.
%
% COMPLEXITY
%   Building A and B costs O(M*N^2) with dense Q; the smallest generalized
%   eigenpair can be obtained efficiently with EIGS for moderate N.
%
% EXAMPLE
%   % Given Q (M-by-N):
%   y = ChooseMonotonicLinearCombination(Q);
%   % Optional: prefer increasing direction
%   D = diff(speye(size(Q,1)));
%   if sum(D*y) < 0, y = -y; end
%
% REFERENCES
%   • Rayleigh–Ritz variational principle for generalized eigenproblems.
%   • Discrete smoothing via difference/Laplacian energies y'*(D'*D)*y.
% See theory in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/LATEX/MonotonicChoiceVectors.pdf
%  JAHO, 23-Oct-2025, Balmes 185, Barcelona
% Comments (and in part implementation) done by ChatGPT 5

 D = diff(speye(size(Q,1)));           % (M-1) x M
A = Q.'*(D.'*D)*Q;
B = Q.'*Q;
[c,~] = eigs(A,B,n,'smallestreal');
y = Q*c;