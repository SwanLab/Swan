classdef LOBPCG
% Testing_LOBPCG
%   Locally Optimized Block Preconditioned Conjugate Gradient implementation
%   for several different examples.

    properties (SetAccess = private)
        K   (:,:) double
        M   (:,:) double
        Kdiag  (:,1) double      % diagonal of K for Jacobi
        Vtrue(:,1) double
        L                        % ichol factor (optional)
        Milu
        Meifem
        Mmult
    end

    properties (SetAccess = public)
        % Solver controls
        b               (1,1) double  = 2         % (# eigenpairs per sweep)
        maxit           (1,1) double  = 20000        % max iterations
        tol             (1,1) double  = 1e-8      % *** relative residual tolerance ***
        use_precond     (1,1) logical = true
        precond_type                 = 'eifem';   % 'ichol' | 'jacobi' | 'none' | 'eifem
        use_physical_x  (1,1) logical = true
        problem         (1,1) double  = 2

        % Print diagnostics
        verbose         (1,1) logical = false
    end

    methods
        function obj = LOBPCG(problem)
            if nargin>0, obj.problem = problem; end
            [obj.K,obj.M] = obj.setup_matrices(obj.problem);

            % Symmetrize once for safety
            obj.K = LOBPCG.symmetrize(obj.K);
            obj.M = LOBPCG.symmetrize(obj.M);

            obj.Kdiag = full(diag(obj.K));

           
        end
    end

    methods
        function results = run_demo(obj)
            
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
            elseif obj.use_precond && strcmpi(obj.precond_type,'eifem')
                    [obj.Mmult, obj.Meifem, obj.Milu] = obj.initEifemPreconditioner();
            end
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
            X = obj.selectInitialGuess(n,b); %select initial subspace
            X = obj.M_orth(X,M);   % ensure X' M X = I -> we do this to keep problem well scaled and numerically stable. 
            % (physically: modes do not share kinetic energy/are
            % dynamically independent. If they are not orthogonal, energy leaks from one mode into another;
            % solver’s subspace becomes coupled and unstable.

            P = [];                % conjugate/search directions
            hist.rnorm  = []; %structure to store residual information
            active = true(1, b);   % all eigenvectors active initially for refinement
            lambda_ritz = zeros(1, b);  % initialize ritz values
            
            % eifem preconditioner initialization
            %LHSfun = @(x) K*x;
            %Milu         = obj.createILUpreconditioner(K);
            %load('eifemPreconditioner.mat')
            %Meifem = @(r) eP.apply(r);
            %Mmult        = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
            
            tic
            for it = 1:obj.maxit
                % 2) Ritz in span(X): best mode approximation in current
                % subspace X
                [lambda_ritz, ~, X] = obj.ritz_step(K, M, X, b);

                % 3) Residual column
                R = K*X(:,active) - M*X(:,active)*diag(lambda_ritz(active));
                
                full_rnorm = obj.computeFullResidual(hist, R, active, b, it,lambda_ritz);                
                hist.rnorm = [hist.rnorm; full_rnorm]; % append residual to structure
  
                newly_converged = (full_rnorm < obj.tol ) & active; % update active set 
                % if any(newly_converged)
                %     fprintf('Locking Mode %s (Converged)\n', ...
                %         mat2str(find(newly_converged)));
                % end 
                active(full_rnorm < obj.tol) = false;

                % Stop if all are converged
                if ~any(active)
                    fprintf('All modes converged.\n');
                    break;
                end

                % 4) Preconditioned residuals Z = P(R), then M-orth proj/out
                % Physical meaning: Nudge each mode in the direction that
                % would annihilate its imbalance 𝐾x = λMx, but only
                % in directions not already spanned by X
                
                Z = obj.apply_prec(R); % preconditioned residuals (steepest-descent corrections)
                Z = obj.M_proj_out(Z, X, M); % remove components of Z already in span(X) so subspace expands
                Z = obj.M_orth(Z, M);

                Xactive = X(:, active); %active modes to refine

                % 5) Build expanded subspace with conjugate directions 
                % P = previous conjugate directions. Gactive defines
                % trial subspace for next ritz projection.
                if isempty(P) 
                    Gactive = [Xactive, Z]; %first iteration no P
                else
                    P = obj.M_proj_out(P, [Xactive Z], M); %remove from P any overlap with [X,Z]
                    P = obj.M_orth(P, M);
                    Gactive = [Xactive, Z, P];
                end
                
                % 6) Ritz on span(Gactive) - expanded subspace; keep best b; update conjugate block
                n_active = sum(active);
                [lam_all, Y, Xactive] = obj.ritz_step(K, M, Gactive, n_active);

                % 7) Partition Y according to Gactive = [Xactive,Z,P]
                hasP = ~isempty(P);
                nX = size(Xactive,2); %size 
                nZ = size(Z,2);
                nP = size(P,2);

                % keep first n_active Ritz vectors (best eigenpairs)
                Ysel = Y(:, 1:n_active);

                %After Ritz step, each new eigenvector X is a combination:
                % Xnew = X·Yx + Z·Yz + P·Yp

                if hasP
                    % Split according to Gactive = [Xactive, Z, P]
                    Yx = Ysel(1:nX, :);
                    Yz = Ysel(nX+1 : nX+nZ, :);
                    Yp = Ysel(nX+nZ+1 : nX+nZ+nP, :);
                    P  = [Z, P] * [Yz; Yp]; %Pnew
                    %in CG new search direction is a combination of new
                    %residual + previous direction. Pnew combines the last
                    %preconditioned residual and old conjugate directions P
                else
                    Yx = Ysel(1:nX, :);
                    Yz = Ysel(nX+1 : nX+nZ, :);
                    P  = Z * Yz;
                end

                %conjugate direction P remembers where the structure was already successfully moving in the last iteration.
                % If we only followed Z each time (pure steepest descent), you would keep “zig-zagging” — overcorrecting the same imbalance back and forth.
                %By forming conjugate directions P, method ensures that the new motion you add does not disturb the energy balance already fixed by previous corrections.

                P = obj.M_proj_out(P, X, M);
                P = obj.M_orth(P, M);

                % 8) Update X and lambdas for active subset
                X(:, active) = Xactive;
                lambda_ritz(active) = lam_all(:).';

                % 9) Soft restart if ill-conditioned (stabilizes convergence)
                Ggram = ([X P].'*(M*[X P]));
                Ggram = LOBPCG.symmetrize(Ggram);
                if rcond(Ggram) < 1e-10 || norm(P,'fro') < 1e-14
                    P = [];
                end
            end
            toc
        end
    end
    %% =======================
    %  Core Solver Stages
    %  =======================

    methods (Access = private)

        function [theta, V, Xout] = ritz_step(obj, K, M, Xin, b)
            % Build an M-orthonormal basis Q for span(Xin) (SVQB fallback inside)
            Q = LOBPCG.M_orth(Xin, M);        % Q' M Q = I -> means every colun has unit generalized mass, and modes are mutually mass-orthogonal
            % Project K; B = I implicitly
            A = Q.' * (K * Q);  A = LOBPCG.symmetrize(A); %modified so more stable numerically than getting eig from A = (X^T K X) & B = (X^T M X)
            
            [V, D] = eig(A, 'vector');
            [theta, p] = sort(real(D), 'ascend');
            V = V(:, p);
            % Lift back and re-orth
            Xout  = Q * V(:, 1:b);
            Xout  = LOBPCG.M_orth(Xout, M);
            theta = theta(1:b).';
        end

        function Z = apply_prec(obj, R, Mmult)
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

                case 'eifem'
                    for i=1:size(R,2)
                        Z(:,i) = obj.Mmult(R(:,i));
                    end
                    
                otherwise
                    Z = R;
            end
        end

        function full_rnorm = computeFullResidual(obj, hist, R, active, b, it,lambda)
            
            rnorm = vecnorm(R, 2, 1);
            if obj.problem == 2
                rnorm =  rnorm ./ max(1, abs(lambda));
            end
            % Store full-length residual history (one row per iteration)
            full_rnorm = zeros(1, b);
            full_rnorm(active) = rnorm;

            % For locked ones, carry previous value (for plotting continuity)
            if it > 1
                prev = hist.rnorm(end, :);
                full_rnorm(~active) = prev(~active);
            else
                full_rnorm(~active) = NaN;  % first iteration
            end
        end

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
                        S = load('PhiCoarseProjectedRed.mat',"PhiFineContRed");
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
    end

    %% =======================
    %  Orthonormalization / Projections (M-inner product)
    %  =======================
    methods (Static, Access = private)
        function X = M_orth(X, M)
            %orthogonality of modes: each mode vibrates independently of
            %others. Mathematically, independence is measured b y how
            %kinetic energies overlap.
            % Mass-weighted orthonormalization: X' M X = I
            if isempty(X), return; end
            G = X.' * (M*X);
            G = LOBPCG.symmetrize(G);
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
            tol = eps(max(s));
            keep = (s > tol);
            U = U(:,keep); s = s(keep);
            if isempty(s)
                % pathological: return zeros to avoid NaNs
                X = zeros(size(X));
                return
            end
            X = X * (U .* (1./sqrt(s)).');
        end

        function Z = M_proj_out(Z, Q, M)
            % M_proj_out removes from new search directions any kinetic-energy overlap with the current subspace,
            % keeps the subspace expanding and the solver converging.
            % Q -> modes we already konw (or subspace they span)
            % Z -> Matrix of new directions to add (proposed corrections oti mprove modes)
            % Goal: remove from Z any component that lies in span(Q), but measured in the M-inner produc
            if isempty(Z) || isempty(Q), return; end
            G = Q.' * (M * Q);  LOBPCG.symmetrize(G); % Compute the Gram matrix and symmetrize
            T = Q.' * (M * Z); % compute cross gram
            Z = Z - Q * (G \ T); % Projecting Z with respect to the mass inner product removes from each correction any overlap in kinetic energy with the existing modes
        end
    end

    %% =======================
    %  Preconditioner Methods
    %  =======================


    methods (Access = private)
        function [Mmult,Meifem, Milu] = initEifemPreconditioner(obj)
            LHSfun = @(x) obj.K*x;
            Milu = obj.createILUpreconditioner(obj.K);
            load('eifemPreconditioner.mat')
            Meifem = @(r) eP.apply(r);
            Mmult = @(r) Preconditioner.multiplePrec(r,LHSfun,Milu,Meifem,Milu);
        end
        function Milu = createILUpreconditioner(obj,LHS)
            s.LHS = LHS;
            s.type = 'ILU';
            M = Preconditioner.create(s);
            Milu = @(r) M.apply(r);
        end
    end

    %% =======================
    %  Problem Setup - Examples
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
                    S1 = load("MassMatrixRed.mat");
                    S2 = load("LHSRed.mat");
                    % Accept either variable names: Mr / M, lhs / LHSr
                    M = S1.Mr;
                    K = S2.LHSr;

                otherwise
                    e = ones(40,1);
                    K = spdiags([-e 2*e -e], -1:1, 4,4);
                    M = spdiags([1; 2; 2; 1.5], 0, 4,4);
            end
        end

    end

    %% =======================
    %  Helpers
    %  =======================

    methods (Static, Access = private)
        function X = symmetrize(X)
            X = (X+X.')/2;
        end
    end
end
