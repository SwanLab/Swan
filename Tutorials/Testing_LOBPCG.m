classdef Testing_LOBPCG
% Testing_LOBPCG
%   Locally Optimized Block Preconditioned Conjugate Gradient implementation
%   for several different examples.

    properties (SetAccess = private)
        K   (:,:) double
        M   (:,:) double
        Kdiag  (:,1) double      % diagonal of K for Jacobi
        Vtrue(:,1) double
        L                        % ichol factor (optional)
    end

    properties (SetAccess = public)
        % Solver controls
        b               (1,1) double  = 2         % (# eigenpairs per sweep)
        maxit           (1,1) double  = 20000        % max iterations
        tol             (1,1) double  = 1e-8      % *** relative residual tolerance ***
        use_precond     (1,1) logical = true
        precond_type                 = 'ichol';   % 'ichol' | 'jacobi' | 'none'
        use_physical_x  (1,1) logical = true
        problem         (1,1) double  = 2

        % Print diagnostics
        verbose         (1,1) logical = false
    end

    methods
        function obj = Testing_LOBPCG(problem)
            if nargin>0, obj.problem = problem; end
            [obj.K,obj.M] = obj.setup_matrices(obj.problem);

            % Symmetrize once for safety
            obj.K = (obj.K+obj.K.')/2;
            obj.M = (obj.M+obj.M.')/2;

            obj.Kdiag = full(diag(obj.K));

            % Try to build an incomplete Cholesky preconditioner if requested
            if obj.use_precond && strcmpi(obj.precond_type,'ichol')
                try
                    opts.type     = 'ict';     % incomplete Cholesky with threshold
                    opts.droptol  = 1e-3;
                    opts.diagcomp = 0.1;       % diagonal compensation
                    obj.L = ichol(obj.K, opts);
                catch
                    warning('ichol failed; falling back to Jacobi preconditioner.');
                    obj.precond_type = 'jacobi';
                    obj.L = [];
                end
            end
        end
    end

    methods
        function results = run_demo(obj)
            % Solve
            [lambda, X, history] = obj.lobpcg_extremal();

            % Package results
            results.lambda   = lambda;
            results.X        = X;
            results.history  = history;

        end
    end

    methods (Access = private)
        function [lambda_ritz, X, hist] = lobpcg_extremal(obj)
            K = obj.K; M = obj.M; b = obj.b;
            n = size(K,1);

            % 1) Initial block
            X = obj.selectInitialGuess(n,b);
            X = obj.M_orth(X,M);   % ensure X' M X = I
            P = [];                % conjugate/search directions
            hist.rnorm  = []; %structure to store residual information
           
            for it = 1:obj.maxit
                % 2) Ritz in span(X): best b modes in current subspace
                [lambda_ritz, ~, X] = obj.ritz_step(K, M, X);

                % 3) Residual column
                R      = K*X - M*X*diag(lambda_ritz);
                rnorm  = vecnorm(R, 2, 1);
                
                hist.rnorm = [hist.rnorm; rnorm];

                if all(rnorm < obj.tol)
                    fprintf('  Converged: all relres < tol.\n');
                    break;
                end

                % 4) Preconditioned residuals Z = P(R), then M-orth proj/out
                Z = obj.apply_prec(R);
                Z = obj.M_proj_out(Z, X, M);
                Z = obj.M_orth(Z, M);

                % 5) Build expanded subspace with conjugate directions
                if isempty(P)
                    G = [X, Z];
                else
                    P = obj.M_proj_out(P, [X Z], M);
                    P = obj.M_orth(P, M);
                    G = [X, Z, P];
                end

                % 6) Ritz on span(G); keep best b; update conjugate block
                [lam_all, Y, X] = obj.ritz_step(K, M, G);

                % Partition Y according to G = [X,Z,P]
                hasP = ~isempty(P);
                Ysel = Y(:,1:b);
                if hasP
                    Yz = Ysel(b+1:2*b, :);
                    Yp = Ysel(2*b+1:3*b, :);
                    P  = [Z, P] * [Yz; Yp];
                else
                    Yz = Ysel(b+1:2*b, :);
                    P  = Z * Yz;
                end
                P = obj.M_proj_out(P, X, M);
                P = obj.M_orth(P, M);

                % Soft restart if ill-conditioned
                Ggram = ([X P].'*(M*[X P])); Ggram = (Ggram+Ggram.')/2;
                if rcond(Ggram) < 1e-10 || norm(P,'fro') < 1e-14
                    P = [];
                end
            end
        end
    end

    methods (Access = private)

        function X = selectInitialGuess(obj,n,b)
            
            if obj.use_physical_x
                if obj.problem == 2
                    % If you have saved Vtrue, you can load and perturb it. Else random.
                    try
                        S = load("Vtrue_prob2.mat","Vtrue");
                        Vtrue = S.Vtrue; %from saved results
                        perturb = 0.9 + 0.1 * (2*rand(size(Vtrue)) - 1);
                        X = perturb(:,1:b) .* Vtrue(:,1:b); %perturb to represent coarse 
                    catch
                        rng(4); X = randn(n,b);
                    end
                elseif obj.problem == 3
                    try
                        S = load('PhiCoarseProjected.mat',"PhiFineCont");
                        X = S.PhiFineCont(:,1:b); %coarse eigenmode projected to fine
                    catch
                        rng(4); X = randn(n,b);
                    end
                else
                    rng(4); X = randn(n,b);
                end
            else
                rng(4); X = randn(n,b);
            end
        end

        function [theta, V, Xout] = ritz_step(obj, K, M, Xin)
            % Build an M-orthonormal basis Q for span(Xin) (SVQB fallback inside)
            Q = Testing_LOBPCG.M_orth(Xin, M);        % Q' M Q = I
            % Project K; B = I implicitly
            A = Q.' * (K * Q);  A = (A+A.')/2;
            % Dense symmetric eig (tiny)
            [V, D] = eig(A, 'vector');
            [theta, p] = sort(real(D), 'ascend');
            V = V(:, p);
            % Lift back and re-orth
            b = obj.b;
            Xout  = Q * V(:, 1:b);
            Xout  = Testing_LOBPCG.M_orth(Xout, M);
            theta = theta(1:b).';
        end

        function Z = apply_prec(obj, R)
            if ~obj.use_precond || strcmpi(obj.precond_type,'none')
                Z = R;
                return
            end
            switch lower(obj.precond_type)
                case 'ichol'
                    if isempty(obj.L)
                        Z = R;  % safety
                    else
                        % Solve (L L^T) Z = R  column-wise
                        Z = obj.L \ (obj.L.' \ R);
                    end
                case 'jacobi'
                    d = obj.Kdiag;
                    d(abs(d) < 1e-14) = 1e-14;
                    Z = R ./ d;      % row-wise scaling (each row / diag(K))
                otherwise
                    Z = R;
            end
        end
    end

    %% =======================
    %  Orthonormalization / projections (M-inner product)
    %  =======================
    methods (Static, Access = private)
        function X = M_orth(X, M)
            % Mass-weighted orthonormalization: X' M X = I
            if isempty(X), return; end
            G = X.' * (M*X);
            G = (G + G.')/2;
            % Try Cholesky first
            [R,p] = chol(G, 'lower');
            if p==0
                X = X / R.';                    % fast path
                return
            end
            % Fallback: SVQB/eigendecomp of Gram
            [U,S] = eig(G);
            s = real(diag(S));
            % keep positive spectrum above tolerance
            tol = max(s)*1e-12;
            keep = (s > tol);
            U = U(:,keep); s = s(keep);
            if isempty(s)
                % pathological: return zeros to avoid NaNs
                X = zeros(size(X));
                return
            end
            X = X * U * diag(1./sqrt(s));
        end

        function Z = M_proj_out(Z, Q, M)
            % Exact M-orthogonal projector even if Q' M Q ≠ I
            if isempty(Z) || isempty(Q), return; end
            G = Q.' * (M * Q);  G = (G+G.')/2;
            T = Q.' * (M * Z);
            Z = Z - Q * (G \ T);
        end
    end

    %% =======================
    %  Problem Setup
    %  =======================
    methods (Access = private)
        function [K,M] = setup_matrices(~,problem)
            switch(problem)
                case 1 % 4-DOF toy
                    e = ones(40,1);
                    K = spdiags([-e 2*e -e], -1:1, 4,4);
                    M = spdiags([1; 2; 2; 1.5], 0, 4,4);

                case 2 % Euler–Bernoulli FE with BCs via DOF elimination
                    ne = 200;  L = 10;  Le = L/ne;
                    nn = ne + 1;  nd = 2*nn;
                    E = 200e9;  I = 1e-6;  rho = 7850;  A = 0.01;

                    Ke = (E*I/Le^3) * [ 12,   6*Le,  -12,   6*Le;
                                         6*Le, 4*Le^2, -6*Le, 2*Le^2;
                                        -12,  -6*Le,   12,  -6*Le;
                                         6*Le, 2*Le^2, -6*Le, 4*Le^2 ];
                    Me = (rho*A*Le/420) * [156,   22*Le,   54,  -13*Le;
                                            22*Le, 4*Le^2, 13*Le, -3*Le^2;
                                            54,    13*Le, 156,   -22*Le;
                                           -13*Le, -3*Le^2, -22*Le, 4*Le^2 ];

                    K = sparse(nd, nd);  M = sparse(nd, nd);
                    for e = 1:ne
                        n1 = e; n2 = e+1;
                        idx = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
                        K(idx,idx) = K(idx,idx) + Ke;
                        M(idx,idx) = M(idx,idx) + Me;
                    end

                    bc_type = "clamped-clamped";
                    switch bc_type
                        case "clamped-clamped", fixed = [1,2, 2*nn-1,2*nn];
                        case "simply-supported", fixed = [1, 2*nn-1];
                        case "clamped-free",     fixed = [1,2];
                        otherwise, error('Unknown bc_type');
                    end
                    free = setdiff(1:(2*nn), fixed);
                    K = K(free, free);
                    M = M(free, free);

                case 3
                    S1 = load("MassMatrix.mat");
                    S2 = load("LHS.mat");
                    % Accept either variable names: Mr / M, lhs / LHSr
                    M = S1.M;
                    K = S2.lhs;

                otherwise
                    e = ones(40,1);
                    K = spdiags([-e 2*e -e], -1:1, 4,4);
                    M = spdiags([1; 2; 2; 1.5], 0, 4,4);
            end
        end

    end
end
