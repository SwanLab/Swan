function [R, t, X_rigid, u_deform] = ExtractRigidMotion_Vectorized(M, X, U)
% Extract the best-fit rigid body motion (rotation and translation)
% from a displacement field using a consistent mass matrix.
%
% INPUTS:
% M  : (3N x 3N) consistent mass matrix
% X  : (N x 3) original coordinates
% U  : (N x 3) displacements
%
% OUTPUTS:
% R        : (3 x 3) rotation matrix
% t        : (1 x 3) translation vector
% X_rigid  : (N x 3) rigid motion coordinates
% u_deform : (N x 3) deformational displacements
% JAHO, 18-May-2025 (+ ChatGPT4)
% See theory in 
% /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN/STRAINING_PART.tex
if nargin == 0
    load('tmp1.mat')
end
ndim = size(X,2) ; 
N = size(X,1);                 % number of nodes
X_def = X + U;                 % deformed configuration

% Build translational matrix T (3N x 3)
T = kron(ones(N,1), eye(ndim));  % 3N x 3

% Flatten coordinates
X_vec = reshape(X',[],1);     % 3N x 1
x_vec = reshape(X_def',[],1); % 3N x 1

% Compute mass-weighted centroids
MX  = T' * M * X_vec;          % 3 x 1
Mx  = T' * M * x_vec;          % 3 x 1
MTT = T' * M * T;              % 3 x 3

X_bar = (MTT \ MX)';           % 1 x 3
x_bar = (MTT \ Mx)';           % 1 x 3

% Centralize configurations
Xc = X - X_bar;                % N x 3
xc = X_def - x_bar;            % N x 3

% Build lumped mass per node
% For consistent M, we compute weights per node via:
% Build lumped mass per node from diagonal blocks of M
 ndof = ndim * N;
w = zeros(N,1);
for d = 1:ndim
    idx = d:ndim:ndof;
    w = w + diag(M(idx, idx));
end
Wmat = spdiags(w, 0, N, N);    % N x N


% Weighted cross-covariance matrix
H = Xc' * Wmat * xc;           % 3 x 3

% SVD to get optimal rotation
[Urot, ~, Vrot] = svd(H);
R = Vrot * Urot';

% Ensure det(R) = +1
if det(R) < 0
    Vrot(:,end) = -Vrot(:,end);
    R = Vrot * Urot';
end

% Optimal translation
t = x_bar - X_bar * R';

% Rigid motion and deformational displacement
X_rigid = X * R' + t;            % N x 3
u_deform = X_def - X_rigid;     % N x 3

end
