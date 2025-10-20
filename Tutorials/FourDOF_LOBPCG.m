classdef FourDOF_LOBPCG
% FOURDOF_LOBPCG
%   Minimal, fully-commented LOBPCG implementation for a 4-DOF spring chain
%   generalized eigenproblem:   K * u = lambda * M * u
%
%   Target: smallest eigenvalues (extremal) using a simple "AMG-like"
%   preconditioner (Jacobi here; plug in your AMG V-cycle in apply_prec).
%
%   Usage:
%       solver = FourDOF_LOBPCG();
%       results = solver.run_demo();
%
%   This class is written to be readable + instructional:
%   - Every method is short and commented
%   - All orthogonality is in the M-inner product
%   - Preconditioner is a single "apply"; swap in AMG where noted
%
%   Author: (you)
%   Date:   (today)

    %% =======================
    %  Properties (data)
    %  =======================
    properties (SetAccess = private)
        % Problem matrices (sparse, SPD)
        K   (:,:) double
        M   (:,:) double

        % Preconditioner data (Jacobi stand-in)
        Kdiag  (:,1) double   % diagonal of K for Jacobi

        % Solver controls
        b        (1,1) double = 2      % block size (# eigenpairs per sweep)
        maxit    (1,1) double = 50     % max outer iterations
        tol      (1,1) double = 1e-10  % residual tolerance (2-norm, per vector)

        % Diagnostics
        verbose  (1,1) logical = true
    end

    %% =======================
    %  Constructor
    %  =======================
    methods
        function obj = FourDOF_LOBPCG()
            % Build the 4-DOF example and a simple preconditioner
            [obj.K,obj.M] = obj.setup_matrices_4dof();
            obj.Kdiag = full(diag(obj.K));
            % --- Safety checks ---
            if isempty(obj.K) || isempty(obj.M)
                error('K or M not initialized properly.');
            end
            if any(obj.Kdiag <= 0)
                error('Jacobi preconditioner invalid: K diagonal must be positive.');
            end
        end
    end

    %% =======================
    %  Public API
    %  =======================
    methods
        function results = run_demo(obj)
            % RUN_DEMO
            %   Complete workflow:
            %     1) Print problem
            %     2) Solve with LOBPCG (preconditioned extremal)
            %     3) Verify vs MATLAB eig (small problem check)
            %     4) Return a struct of results

            if obj.verbose
                obj.print_header();
                obj.print_problem();
            end

            % Solve
            [lambda, X, history] = obj.lobpcg_extremal();

            % Truth via dense eig (for demonstration/verification only)
            [Vtrue, Dtrue] = eig(full(obj.K), full(obj.M));
            [lam_true, p]  = sort(diag(Dtrue), 'ascend');
            Vtrue = Vtrue(:, p);

            if obj.verbose
                fprintf('\n=== Converged eigenvalues (LOBPCG, smallest %d) ===\n', obj.b);
                disp(lambda.');

                fprintf('--- Reference (eig, smallest %d) ---\n', obj.b);
                disp(lam_true(1:obj.b).');

                % Residual checks
                rnorms = vecnorm(obj.K*X - obj.M*X.*lambda.', 2, 1);
                fprintf('Residual 2-norms per eigenpair:\n');
                disp(rnorms);

                % M-orthogonality check
                MX = obj.M*X;
                Gram = X.'*MX;
                fprintf('X^T M X (should be ~I):\n');
                disp(Gram);
            end

            % Package results
            results.lambda   = lambda;
            results.X        = X;
            results.history  = history;
            results.lambda_ref = lam_true;
            results.V_ref      = Vtrue;
        end
    end

    %% =======================
    %  Core solver (LOBPCG)
    %  =======================
    methods (Access = private)
        function [lambda, X, hist] = lobpcg_extremal(obj)
            % LOBPCG_EXTREMAL
            %   Block LOBPCG for smallest eigenvalues of K u = lambda M u.
            %   Preconditioner P ~ K^{-1} is applied as a black-box "apply"
            %   to the block residuals (Jacobi here, replace with AMG).
            %
            %   Returns:
            %     lambda  : 1xb vector of Ritz values (ascending)
            %     X       : nxb M-orthonormal Ritz vectors
            %     hist    : struct with iteration diagnostics

            K = obj.K; M = obj.M; b = obj.b;
            n = size(K,1);
            if isempty(K) || isempty(M)
                error('K and M must be initialized before running LOBPCG.');
            end
            if isempty(obj.Kdiag)
                obj.Kdiag = full(diag(K)); % safety net
            end

            % --- Initialization: random block, then M-orthonormalize
            rng(4);
            X = randn(n, b);
            X = obj.M_orth(X, M);
            W = [];  % conjugate block (None on first sweep)

            hist.rnorms = [];

            if obj.verbose
                fprintf('\n>> LOBPCG (extremal) start: b=%d, tol=%.1e, maxit=%d\n', ...
                    b, obj.tol, obj.maxit);
            end

            for it = 1:obj.maxit
                % 1) Ritz on current block (small projected GEP)
                [lambda, X] = obj.ritz_step(K, M, X);

                % 2) Residuals
                R = K*X - M*X.*lambda.';   % columns are residuals r_i
                rnorm = vecnorm(R, 2, 1);
                hist.rnorms = [hist.rnorms; rnorm]; %#ok<AGROW>

                if obj.verbose
                    fprintf('  it=%2d | lambda=', it);
                    fprintf(' %.6f', lambda);
                    fprintf(' | ||r||2=');
                    fprintf(' %.2e', rnorm);
                    fprintf('\n');
                end

                % Convergence check
                if all(rnorm < obj.tol)
                    if obj.verbose
                        fprintf('  Converged: all residuals < tol.\n');
                    end
                    break;
                end

                % 3) Precondition residuals: Z = P(R)
                Z = obj.apply_prec(R);     % <-- AMG V-cycle goes here
                Z = obj.M_proj_out(Z, X, M);
                Z = obj.M_orth(Z, M);

                % 4) Conjugate directions block W (optional but standard)
                if ~isempty(W)
                    W = obj.M_proj_out(W, [X Z], M);
                    W = obj.M_orth(W, M);
                    G = [X, Z, W];
                else
                    G = [X, Z];
                end

                % 5) Rayleigh–Ritz on expanded subspace and keep best b
                AG = G.'*(K*G);
                BG = G.'*(M*G);
                [V,D] = eig(AG, BG);
                [lam_all, idx] = sort(diag(D), 'ascend');
                V = V(:, idx);

                X = G * V(:,1:b);
                X = obj.M_orth(X, M);
                lambda = lam_all(1:b).';

                % Next W: next set of Ritz vectors from G (or just Z)
                if size(V,2) >= 2*b
                    W = G * V(:, b+1:2*b);
                else
                    W = Z;
                end
            end
        end
    end

    %% =======================
    %  Linear algebra helpers
    %  =======================
    methods (Access = private)
        function [theta, Xout] = ritz_step(~, K, M, Xin)
            % RITZ_STEP
            %   Small projected generalized eigenproblem on span(Xin).
            A = Xin.' * (K*Xin);
            B = Xin.' * (M*Xin);
            [V,D] = eig(A,B);
            [theta, p] = sort(diag(D), 'ascend');
            Xout = Xin * V(:,p);
            Xout = FourDOF_LOBPCG.M_orth(Xout, M); % static ok
            theta = theta.';
        end

        function Z = apply_prec(obj, R)
            % APPLY_PREC
            %   Preconditioner apply: Z = P(R)
            %
            %   Here: Jacobi (diagonal) preconditioner as a placeholder for
            %   an AMG V-cycle. Replace the body with a call to your AMG.
            %
            %   Example replacement:
            %       Z = AMG_Vcycle_apply(obj.AMG, R);
            %
            %   Note: For a robust AMG on elasticity, include RBMs as
            %   near-nullspace vectors in SA-AMG setup.
            Z = R ./ obj.Kdiag;   % element-wise divide each row by diag(K)
        end
    end

    %% =======================
    %  Orthonormalization / projections (M-inner product)
    %  =======================
    methods (Static, Access = private)
        function X = M_orth(X, M)
            % M_ORTH
            %   Mass-weighted orthonormalization: X' M X = I
            G = X.' * (M*X);
            % Cholesky safer: assert SPD
            R = chol((G+G.')/2, 'lower');  % symmetrize for numerical safety
            X = X / R.';
        end

        function Z = M_proj_out(Z, Q, M)
            % M_PROJ_OUT
            %   Project Z to be M-orthogonal to the columns of Q
            if isempty(Q), return; end
            Z = Z - Q * (Q.' * (M*Z));
        end
    end

    %% =======================
    %  Problem setup & printing
    %  =======================
    methods (Access = private)
        function [K,M] = setup_matrices_4dof(obj)
            % SETUP_MATRICES_4DOF
            %   4-DOF spring chain, fixed at the right end:
            %     K = tridiag([ -1 2 -1 ]), M diagonal with varying masses
            e = ones(4,1);
            K = spdiags([-e 2*e -e], -1:1, 4,4);
            M = spdiags([1; 2; 2; 1.5], 0, 4,4);


        end

        function print_header(~)
            fprintf('=====================================================\n');
            fprintf('  LOBPCG (preconditioned extremal) on 4-DOF example \n');
            fprintf('  Problem: K u = lambda M u (SPD)\n');
            fprintf('=====================================================\n');
        end

        function print_problem(obj)
            fprintf('\n-- Matrices (sparse display) --\n');
            disp(obj.K);
            disp(obj.M);

            % Reference: true smallest few eigenvalues
            [~, D] = eig(full(obj.K), full(obj.M));
            lam = sort(diag(D), 'ascend');
            fprintf('Reference eigenvalues (all):\n');
            disp(lam.');
        end
    end
end
