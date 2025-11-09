classdef Testing_LOBPCG    
% Testing_LOBPCG
%   Locally Optimized Block Preconditioned Conjugate Gradient implementation
%   for several different examples.

    properties (SetAccess = private)

        K   (:,:) double
        M   (:,:) double
        Kdiag  (:,1) double   % diagonal of K for Jacobi
        Vtrue(:,1) double
        
    end

    properties (SetAccess = public)
        % Solver controls
        b             (1,1) double = 2      % (# eigenpairs per sweep)
        maxit         (1,1) double = 50     % max iterations
        tol           (1,1) double = 1e-10  % converg criteria
        use_precond   (1,1) logical = true  % use preconditioner (true) or simply I for residual (false)
        use_physical_x (1,1) logical = false %use physics derived initial aproximation for ritz step

        % Print diagnostics
        verbose  (1,1) logical = true

    end

    methods
        function obj = Testing_LOBPCG(problem)

            [obj.K,obj.M] = obj.setup_matrices(problem);
            obj.Kdiag = full(diag(obj.K)); %simple preconditioner
            
        end
    end


    methods
        function results = run_demo(obj)

            if obj.verbose
                obj.print_header();
                obj.print_problem();
            end

             % true values for comparison
            [Vtrue, Dtrue] = eigs(full(obj.K), full(obj.M), obj.b, "smallestabs");
            [lam_true, p]  = sort(diag(Dtrue), 'ascend');
            Vtrue = Vtrue(:, p);
            save('Vtrue_prob2.mat','Vtrue')
            save('lam_true_prob2.mat','lam_true')
            % Solve
            [lambda, X, history] = obj.lobpcg_extremal();

            % Package results
            results.lambda   = lambda;
            results.X        = X;
            results.history  = history;
            results.lambda_ref = lam_true;
            results.V_ref      = Vtrue;
        end
    end

    methods (Access = private)
        function [lambda_ritz, X, hist] = lobpcg_extremal(obj) %  Core solver (LOBPCG)
            % LOBPCG_EXTREMAL
            %   LOBPCG for smallest eigenvalues of K u = lambda M u.
            %   Preconditioner P ~ K^{-1} is applied as a black-box "apply"
            %   to the block residuals (Jacobi here, replace with AMG).
            %
            %   Returns:
            %     lambda  : 1xb vector of Ritz values (ascending)
            %     X       : nxb M-orthonormal Ritz vectors
            %     hist    : struct with iteration diagnostics

            K = obj.K; M = obj.M; b = obj.b;
            n = size(K,1);
            
            % --- Initialization: random block, then M-orthonormalize
            if obj.use_physical_x == true
                load("Vtrue_prob2.mat");
                %perturbation = 0.9;
                perturbation = 0.9 + 0.1 * (2*rand(size(Vtrue)) - 1);  % uniform in [0.8, 1.0]
                X = perturbation .* Vtrue;

            else
                rng(4);
                X = randn(n, b); %we use reduced basis X random here but ideally this would be derived from physics and be made up of columns of coarse eigenvectors.
                %X = obj.M_orth(X, M);
            end
            X = obj.M_orth(X,M);
            P = [];  % <-- previous search directions (the conjugate block)

            hist.rnorms = [];

            tic
            for it = 1:obj.maxit
                % 1) Ritz on current block (small projected GEP)
                [lambda_ritz,u_ritz, X] = obj.ritz_step(K, M, X);

                % 2) Residuals
                R = K*X - M*X*diag(lambda_ritz);  
                rnorm = vecnorm(R, 2, 1);
                relres = rnorm ./ ( vecnorm(K*X,2,1) + abs(lambda_ritz).*vecnorm(M*X,2,1) );
                hist.rnorms = [hist.rnorms; relres];

                if all(relres < 1e-8)
                    fprintf('  Converged: all residuals < rel tol.\n');
                    break;
                end

                % 3) Precondition residuals: Z = P(R)
                Z = obj.apply_prec(R);     % <-- preconditioning cycle
                Z = obj.M_proj_out(Z, X, M);
                Z = obj.M_orth(Z, M);

                % 4) Conjugate directions block P (previous search directions)
                %    Build expanded subspace G = [X, Z, P]
                if isempty(P)
                    G = [X, Z];
                else
                    % keep P independent of current [X,Z]
                    P = obj.M_proj_out(P, [X Z], M);
                    P = obj.M_orth(P, M);
                    G = [X, Z, P];
                end

                % 5) Rayleigh–Ritz on expanded subspace and keep best b
                %    ritz_step returns: (theta_all, Y, X_new) with Y the eigenvectors in basis G
                [lam_all, Y, X] = obj.ritz_step(K, M, G);  % X is already M-orthonormal

                % Partition Y consistently with G = [X, Z, P]
                % NOTE: sizes adapt if P is empty
                bx = size(X,2);  % but X here is already updated; use block sizes from G:
                bX = size(G,2);  % total columns in G
                b  = obj.b;
                hasP = (~isempty(P));

                if hasP
                    % sizes: [X,Z,P] = [b, b, b] in your current setup
                    Ysel = Y(:,1:b);          % coefficients of selected Ritz vectors
                    Yx = Ysel(1:b,         :);
                    Yz = Ysel(b+1:2*b,     :);
                    Yp = Ysel(2*b+1:3*b,   :);
                    % Form new conjugate block from components that live in span{Z,P}
                    P = [Z, P] * [Yz; Yp];
                else
                    Ysel = Y(:,1:b);
                    Yx = Ysel(1:b,     :);
                    Yz = Ysel(b+1:2*b, :);
                    % First iteration: no P yet; just carry Z-part
                    P = Z * Yz;
                end

                % Make P M-orthogonal to the new X, and well-conditioned
                P = obj.M_proj_out(P, X, M);
                P = obj.M_orth(P, M);

                % (Optional) Soft-restart if [X,P] Gram is ill-conditioned
                Ggram = ( [X, P].' * (M*[X, P]) ); Ggram = (Ggram+Ggram.')/2;
                if rcond(Ggram) < 1e-10 || norm(P,'fro') < 1e-14
                    P = [];  % drop directions to avoid stagnation/instability
                end
            end
            toc
        end
    end

    %% =======================
    %  Linear algebra helpers
    %  =======================
    methods (Access = private)
        function [theta, V, Xout] = ritz_step(obj, K, M, Xin)
            % RITZ_STEP
            %   Small projected generalized eigenproblem on span(Xin).
            A = Xin.' * (K*Xin); %reduced "condensed" stiffness matrix\equa
            B = Xin.' * (M*Xin); %reduced "condensed" mass matrix
            A = (A+A.')/2; B = (B+B.')/2;        % symmetrize for safety
            % A y = lambda B y (ritz problem)
            [V,D] = eig(A,B,'vector');
            [theta, p] = sort(real(D), 'ascend');
            V = V(:,p);
            Xout = Xin * V(:,1:obj.b);
            Xout = Testing_LOBPCG.M_orth(Xout, M); % static ok
            theta = theta(1:obj.b).';
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
            if obj.use_precond == true
                d = obj.Kdiag;
                d(abs(d) < 1e-14) = 1e-14;
                Z = R ./ d; % element-wise divide each row by diag(K)
            else
                Z = R;
            end
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

    methods (Access = private)
        function [K,M] = setup_matrices(obj,problem)

            switch(problem)
                case 1 % SETUP_MATRICES_4DOF- 4-DOF spring chain, fixed at the right end:
                    e = ones(40,1);
                    K = spdiags([-e 2*e -e], -1:1, 4,4);
                    M = spdiags([1; 2; 2; 1.5], 0, 4,4);
                case 2 % Euler–Bernoulli FE with BCs handled by DOF elimination
                    % Discretization
                    ne = 200;         % number of elements (adjust as you like)
                    L  = 10;
                    Le = L/ne;
                    nn = ne + 1;      % number of nodes
                    nd = 2*nn;        % DOFs: [w1,theta1, w2,theta2, ..., w_nn,theta_nn]

                    % Beam parameters
                    E = 200e9;  I = 1e-6;
                    rho = 7850; A = 0.01;

                    % Element stiffness & mass (Hermite cubic)
                    Ke = (E*I/Le^3) * [ 12,   6*Le,  -12,   6*Le;
                        6*Le, 4*Le^2, -6*Le, 2*Le^2;
                        -12,  -6*Le,   12,  -6*Le;
                        6*Le, 2*Le^2, -6*Le, 4*Le^2 ];

                    Me = (rho*A*Le/420) * [156,   22*Le,   54,  -13*Le;
                        22*Le, 4*Le^2, 13*Le, -3*Le^2;
                        54,    13*Le, 156,   -22*Le;
                        -13*Le, -3*Le^2, -22*Le, 4*Le^2 ];

                    % Global assembly
                    K = sparse(nd, nd);
                    M = sparse(nd, nd);
                    for e = 1:ne
                        % dof map for element e
                        n1 = e; n2 = e+1;
                        idx = [ 2*n1-1, 2*n1, 2*n2-1, 2*n2 ];
                        K(idx,idx) = K(idx,idx) + Ke;
                        M(idx,idx) = M(idx,idx) + Me;
                    end

                    % === Pick BC type ===
                    % "clamped-clamped": w=0 and theta=0 at both ends
                    % "simply-supported": w=0 at both ends (rotations free)
                    % "clamped-free": cantilever (left clamped, right free)
                    bc_type = "clamped-clamped";

                    switch bc_type
                        case "clamped-clamped"
                            fixed = [1, 2,  2*nn-1, 2*nn];     % [w1,theta1, wN,thetaN]

                        case "simply-supported"
                            fixed = [1, 2*nn-1];               % constrain only deflection at ends

                        case "clamped-free"                    % cantilever
                            fixed = [1, 2];                    % clamp left end; right end free

                        otherwise
                            error('Unknown bc_type');
                    end

                    % Eliminate constrained DOFs
                    free = setdiff(1:nd, fixed);
                    K = K(free, free);
                    M = M(free, free);

                    % (Optional) Symmetrize for numerical safety
                    K = (K+K.')/2;  M = (M+M.')/2;

                case 3
                    load("MassMatrixRed.mat")
                    load("LHSRed.mat")
                    K = LHSr;
                    M = Mr;
                case default
                    e = ones(40,1);
                    K = spdiags([-e 2*e -e], -1:1, 4,4);
                    M = spdiags([1; 2; 2; 1.5], 0, 4,4);
            end
        end

        function print_header(~)
            fprintf('=====================================================\n');
            fprintf('  LOBPCG on 4-DOF example \n');
            fprintf('  Problem: K u = lambda M u (SPD)\n');
            fprintf('=====================================================\n');
        end

        function print_problem(obj)
            fprintf('\n-- K & M Matrices ) --\n');
            disp(full(obj.K));
            disp(full(obj.M));

            [~, D] = eigs(full(obj.K), full(obj.M), obj.b, "smallestabs");
            lam = sort(diag(D), 'ascend');
            fprintf('Reference eigenvalues (all):\n');
            disp(lam.');
        end
    end
end
