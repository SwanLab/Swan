classdef TutorialMixedStokes2D < handle
    % Stokes 2D on the unit square, mixed formulation P1-P1 (triangles)
    %
    %   -nu Lap(u) + grad(p) = f,   div(u) = 0   in (0,1)^2
    %    u = 0 on the boundary,     p(0,0) = -1/2
    %
    % Weak form (pressure gradient kept on the trial side):
    %   nu (grad u, grad du) + (grad p, du) = (f, du)
    %   (grad dp, u)                        = 0       (= -(dp, div u))
    %
    % Manufactured solution (div u = 0, u = 0 on the boundary):
    %   psi = x^2 (1-x)^2 y^2 (1-y)^2,  u = (dpsi/dy, -dpsi/dx)
    %   p   = x^3 + y^3 - 1/2
    %
    % Cases (same as the 1D tutorial):
    %   MixedFormulation_classic        Galerkin, no stabilization (singular)
    %   MixedFormulation_bubble         P1 + cubic bubble / P1 (MINI), exact condensation
    %   MixedFormulation_eigen(type)    fine space {b*l1, b*l2, b*l3} per component,
    %                                   K_vv ~ tau^-1 M_vv, tau from eig(K_vv, M_vv)
    %                                   type = 'min' | 'mean' | 'max' | 'exact'

    properties (Access = private)
        mesh
        nu
        uFun, pFun
        bcU, bcP
        tau
    end

    methods (Access = public)

        function obj = TutorialMixedStokes2D()
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()
            obj.createBoundaryConditions()

            obj.MixedFormulation_classic()
            obj.computeNormError('u-p (P1-P1) no stab.');
            obj.plotResults('u-p (P1-P1) no stab.');

            obj.MixedFormulation_bubble()
            obj.computeNormError('u-p (P1-P1) bubble funct.');
            obj.plotResults('u-p (P1-P1) bubble funct.');
            obj.plotSubscale('bubble');
            uB = double(obj.uFun.fValues);  pB = double(obj.pFun.fValues);

            types = {'exact', 'min', 'mean', 'max'};
            for i = 1:numel(types)
                obj.MixedFormulation_eigen(types{i})
                obj.computeNormError(['u-p (P1-P1) tau ', types{i}]);
            end

            obj.MixedFormulation_eigen('min')
            uE = double(obj.uFun.fValues);  pE = double(obj.pFun.fValues);
            obj.plotResults('u-p (P1-P1) tau min');
            obj.plotSubscale('min');

            fprintf('max|u_bubble - u_tauMin| = %e\n', max(abs(uB(:) - uE(:))));
            fprintf('max|p_bubble - p_tauMin| = %e\n', max(abs(pB(:) - pE(:))));

            obj.plotKvvComparison(1);         % element 1: prints K_vv, tau^-1 M_vv and S_pp
            % obj.plotKvvComparison(1, true); % same, plus the figure
        end

    end

    methods (Access = private)

        function createMesh(obj)
            n = 10;                                  % n x n nodes (+ one node per cell centre)
            obj.mesh = UnitTriangleMesh(n, n);
        end

        function createMaterialProperties(obj)
            viscosity = 1;
            obj.nu = ConstantFunction.create(viscosity, obj.mesh);
        end

        function createSolutionField(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh, 2, 'P1');
            obj.pFun = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end

        function createBoundaryConditions(obj)
            obj.createBoundaryConditionsVelocity();
            obj.createBoundaryConditionsPressure();
        end

        function createBoundaryConditionsVelocity(obj)
            TOL = 1e-12;
            onBoundary = @(coor) abs(coor(:,1)) < TOL | abs(coor(:,1)-1) < TOL | ...
                                 abs(coor(:,2)) < TOL | abs(coor(:,2)-1) < TOL;

            % Velocity: u = 0 on the whole boundary, both components
            sDir{1}.domain    = onBoundary;
            sDir{1}.direction = [1 2];
            sDir{1}.value     = 0;
            sDir{1}.ndim      = length(sDir{1}.direction);

            obj.bcU = obj.createBoundaryConditionsFromInput(sDir);
        end

        function createBoundaryConditionsPressure(obj)
            TOL = 1e-12;
            isCorner = @(coor) abs(coor(:,1)) < TOL & abs(coor(:,2)) < TOL;

            % Pressure: p(0,0) = -1/2 (pressure defined up to a constant)
            sDir{1}.domain    = isCorner;
            sDir{1}.direction = [1];
            sDir{1}.value     = -0.5;
            sDir{1}.ndim      = length(sDir{1}.direction);

            obj.bcP = obj.createBoundaryConditionsFromInput(sDir);
        end

        function bc = createBoundaryConditionsFromInput(obj, sDir)
            dirichletFun = [];
            for i = 1:numel(sDir)
                dir = DirichletCondition(obj.mesh, sDir{i});
                dirichletFun = [dirichletFun, dir];
            end
            s.dirichletFun = dirichletFun;
            s.periodicFun  = [];
            s.pointloadFun = [];
            s.mesh         = obj.mesh;
            bc = BoundaryConditions(s);
        end

        % ------------------------------------------------------------------
        %  Cases
        % ------------------------------------------------------------------

        function MixedFormulation_classic(obj)
            LHS = obj.createLHS_v1();
            RHS = obj.createRHS_v1();
            [LHSr, RHSr] = obj.reduceSystem(LHS, RHS);
            fprintf('Galerkin P1-P1: rcond = %.2e\n', rcond(full(LHSr)));
            x = pinv(full(LHSr))*RHSr;               % singular: minimum-norm solution
            obj.setSolution(x);
        end

        function MixedFormulation_bubble(obj)
            LHS = obj.createLHS_v1();
            RHS = obj.createRHS_v1();
            [Spp, Fp] = obj.createSubscaleCondensation('bubble');
            [LHS, RHS] = obj.addPressureBlock(LHS, RHS, Spp, Fp);
            [LHSr, RHSr] = obj.reduceSystem(LHS, RHS);
            obj.setSolution(LHSr\RHSr);
        end

        function MixedFormulation_eigen(obj, tauType)
            LHS = obj.createLHS_v1();
            RHS = obj.createRHS_v1();
            [Spp, Fp] = obj.createSubscaleCondensation(tauType);
            [LHS, RHS] = obj.addPressureBlock(LHS, RHS, Spp, Fp);
            [LHSr, RHSr] = obj.reduceSystem(LHS, RHS);
            obj.setSolution(LHSr\RHSr);
        end

        % ------------------------------------------------------------------
        %  Galerkin system
        % ------------------------------------------------------------------

        function LHS = createLHS_v1(obj)
            nu = obj.nu;
            % nu (grad u, grad du)
            K = IntegrateLHS(@(u,v) nu.*DDP(Grad(v),Grad(u)), obj.uFun, obj.uFun, obj.mesh, 'Domain', 2);
            % (grad p, du)   (= -(p, div du) since du = 0 on the boundary)
            B = IntegrateLHS(@(p,v) DP(v,Grad(p)), obj.uFun, obj.pFun, obj.mesh, 'Domain', 2);
            Z = sparse(obj.pFun.nDofs, obj.pFun.nDofs);
            LHS = [K, B; B', Z];
        end

        function RHS = createRHS_v1(obj)
            f  = AnalyticalFunction.create(@(X) obj.loadVector(X), obj.mesh);
            Fu = IntegrateRHS(@(v) DP(v,f), obj.uFun, obj.mesh, 'Domain', 4);
            Fp = zeros(obj.pFun.nDofs, 1);
            RHS = [Fu; Fp];
        end

        function fX = loadVector(obj, X)
            x = X(1,:,:);  y = X(2,:,:);
            [fx, fy] = obj.load(x, y);
            fX = [fx; fy];
        end

        function [fx, fy] = load(obj, x, y)
            nu = obj.nu.constant;
            fx = -24*nu*x.^4.*y + 12*nu*x.^4 + 48*nu*x.^3.*y - 24*nu*x.^3 - 48*nu*x.^2.*y.^3 ...
                 + 72*nu*x.^2.*y.^2 - 48*nu*x.^2.*y + 12*nu*x.^2 + 48*nu*x.*y.^3 - 72*nu*x.*y.^2 ...
                 + 24*nu*x.*y - 8*nu*y.^3 + 12*nu*y.^2 - 4*nu*y + 3*x.^2;
            fy = 48*nu*x.^3.*y.^2 - 48*nu*x.^3.*y + 8*nu*x.^3 - 72*nu*x.^2.*y.^2 + 72*nu*x.^2.*y ...
                 - 12*nu*x.^2 + 24*nu*x.*y.^4 - 48*nu*x.*y.^3 + 48*nu*x.*y.^2 - 24*nu*x.*y + 4*nu*x ...
                 - 12*nu*y.^4 + 24*nu*y.^3 - 12*nu*y.^2 + 3*y.^2;
        end

        % ------------------------------------------------------------------
        %  Subscales: bubble (exact condensation) or tau from eigenvalues
        % ------------------------------------------------------------------

        function [Spp, Fp] = createSubscaleCondensation(obj, type)
            % Per element, fine velocity modes phi_k e_c (c = x, y):
            %   'bubble'      : phi = 27 l1 l2 l3                  (1 mode/comp.)
            %   eigen types   : phi = b l1, b l2, b l3, b = l1 l2 l3 (3 modes/comp.)
            % Blocks (velocity test rows, as in createLHS_v1):
            %   Kvv = nu (grad phi_k, grad phi_l)   (block diagonal in c)
            %   Bvp = (dN_j/dx_c, phi_k)            (grad p, dv)
            %   fv  = (f_c, phi_k)
            %   Kuv = nu (grad phi, grad N) = 0     (P1 gradient constant, phi = 0 on dK)
            % Condensation:  Spp = -Bvp' Kvv^-1 Bvp,  Fp = -Bvp' Kvv^-1 fv
            % tau version:   Kvv ~ tau^-1 Mvv
            nu     = obj.nu.constant;
            coord  = obj.mesh.coord;
            connec = obj.mesh.connec;
            nElem  = size(connec, 1);
            nP     = obj.pFun.nDofs;
            [L, w] = obj.triangleQuadrature(6);     % barycentric coords (nG x 3), weights (sum 1/2)

            I = zeros(9*nElem,1); J = I; V = I;
            Fp = zeros(nP, 1);
            obj.tau = zeros(nElem, 1);
            for e = 1:nElem
                nodes = connec(e,:);
                xe = coord(nodes,:);                         % 3 x 2
                T  = [ones(3,1), xe];
                C  = T \ eye(3);                             % columns: coefficients of l_i
                dL = C(2:3,:);                               % 2 x 3, grad l_i
                A2 = abs(det(T));                            % 2 * area
                xg = L*xe;                                   % Gauss points (nG x 2)
                [fx, fy] = obj.load(xg(:,1), xg(:,2));

                [phi, dphi] = obj.fineModes(L, dL, type);    % nG x m, nG x m x 2
                m  = size(phi, 2);
                W  = A2*w;                                   % weights x detJ

                K1 = zeros(m); M1 = zeros(m);
                for d = 1:2
                    K1 = K1 + nu*dphi(:,:,d)'*(W.*dphi(:,:,d));
                end
                M1 = phi'*(W.*phi);
                Kvv = blkdiag(K1, K1);                       % [x-modes, y-modes]
                Mvv = blkdiag(M1, M1);
                Bvp = [ (phi'*W)*dL(1,:) ;                    % x-modes: dN_j/dx * int phi
                        (phi'*W)*dL(2,:) ];                   % y-modes: dN_j/dy * int phi
                fv  = [ phi'*(W.*fx) ; phi'*(W.*fy) ];

                switch type
                    case {'bubble', 'exact'}
                        Se = -Bvp'*(Kvv\Bvp);
                        Fe = -Bvp'*(Kvv\fv);
                    otherwise
                        lam = eig(K1, M1);
                        switch type
                            case 'max';  t = 1/max(lam);
                            case 'mean'; t = 2/(min(lam) + max(lam));
                            case 'min';  t = 1/min(lam);
                        end
                        Se = -t*Bvp'*(Mvv\Bvp);
                        Fe = -t*Bvp'*(Mvv\fv);
                        obj.tau(e) = t;
                end

                k = 9*(e-1) + (1:9);
                [jj, ii] = meshgrid(nodes, nodes);
                I(k) = ii(:); J(k) = jj(:); V(k) = Se(:);
                Fp(nodes) = Fp(nodes) + Fe;
            end
            Spp = sparse(I, J, V, nP, nP);
        end

        function [phi, dphi] = fineModes(~, L, dL, type)
            % L: nG x 3 barycentric coordinates, dL: 2 x 3 gradients of l_i
            l1 = L(:,1); l2 = L(:,2); l3 = L(:,3);
            b  = l1.*l2.*l3;
            db = zeros(numel(l1), 2);
            for d = 1:2
                db(:,d) = dL(d,1)*l2.*l3 + dL(d,2)*l1.*l3 + dL(d,3)*l1.*l2;
            end
            if strcmp(type, 'bubble')
                phi  = 27*b;
                dphi = reshape(27*db, [], 1, 2);
            else
                phi  = [b.*l1, b.*l2, b.*l3];
                dphi = zeros(numel(l1), 3, 2);
                for i = 1:3
                    for d = 1:2
                        dphi(:,i,d) = db(:,d).*L(:,i) + b*dL(d,i);
                    end
                end
            end
        end

        function [L, w] = triangleQuadrature(obj, n)
            % Collapsed (Duffy) Gauss rule on the reference triangle,
            % exact for polynomials up to degree 2n-2. Weights sum to 1/2.
            [t, wt] = obj.gaussLegendre01(n);
            [U, Vv] = meshgrid(t, t);
            [WU, WV] = meshgrid(wt, wt);
            xi  = U(:);
            eta = Vv(:).*(1 - xi);
            w   = WU(:).*WV(:).*(1 - xi);
            L   = [1 - xi - eta, xi, eta];
        end

        function [x, w] = gaussLegendre01(~, m)
            % Gauss-Legendre points and weights on [0,1] (Golub-Welsch)
            k  = 1:m-1;
            be = k./sqrt(4*k.^2 - 1);
            [Vg, Dg] = eig(diag(be,1) + diag(be,-1));
            [x, idx] = sort((diag(Dg) + 1)/2);
            w  = Vg(1,idx)'.^2;                       % sums to 1 on [0,1]
        end

        % ------------------------------------------------------------------
        %  Assembly helpers
        % ------------------------------------------------------------------

        function [LHS, RHS] = addPressureBlock(obj, LHS, RHS, Spp, Fp)
            nUt = obj.uFun.nDofs;
            iP  = nUt + (1:obj.pFun.nDofs);
            LHS(iP,iP) = LHS(iP,iP) + Spp;
            RHS(iP)    = RHS(iP)    + Fp;
        end

        function [LHSr, RHSr] = reduceSystem(obj, LHS, RHS)
            nUt   = obj.uFun.nDofs;
            dofFU = obj.bcU.free_dofs(:);
            dofFP = obj.bcP.free_dofs(:);
            free  = [dofFU; nUt + dofFP];
            dir   = [obj.bcU.dirichlet_dofs(:); nUt + obj.bcP.dirichlet_dofs(:)];
            vals  = [obj.bcU.dirichlet_vals(:); obj.bcP.dirichlet_vals(:)];
            LHSr  = LHS(free, free);
            RHSr  = RHS(free) - LHS(free, dir)*vals;
        end

        function setSolution(obj, x)
            dofFU = obj.bcU.free_dofs(:);
            dofFP = obj.bcP.free_dofs(:);
            nU = length(dofFU);
            nP = length(dofFP);

            uT = zeros(obj.uFun.nDofs,1);
            uT(obj.bcU.dirichlet_dofs) = obj.bcU.dirichlet_vals;
            uT(dofFU,1) = x(1:nU);
            uT = reshape(uT, obj.uFun.ndimf, [])';
            obj.uFun.setFValues(uT);

            pT = zeros(obj.pFun.nDofs,1);
            pT(obj.bcP.dirichlet_dofs) = obj.bcP.dirichlet_vals;
            pT(dofFP,1) = x(nU+1:nU+nP);
            obj.pFun.setFValues(pT);
        end

        % ------------------------------------------------------------------
        %  Errors and plots
        % ------------------------------------------------------------------

        function [errUH1, errPL2] = computeNormError(obj, label)
            % |u - u_h|_H1 and ||p - p_h||_L2, integrated element by element
            coord  = obj.mesh.coord;
            connec = obj.mesh.connec;
            uh = double(obj.uFun.fValues);                   % nNodes x 2
            ph = double(obj.pFun.fValues);
            [L, w] = obj.triangleQuadrature(5);
            eU = 0; eP = 0;
            for e = 1:size(connec,1)
                nodes = connec(e,:);
                xe = coord(nodes,:);
                T  = [ones(3,1), xe];
                C  = T \ eye(3);
                dL = C(2:3,:);
                W  = abs(det(T))*w;
                xg = L*xe;  x = xg(:,1);  y = xg(:,2);

                Guh = dL*uh(nodes,:);                        % [dux/dx duy/dx; dux/dy duy/dy]
                duxdx = 4.*x.*y.*(x - 1).*(2.*x - 1).*(y - 1).*(2.*y - 1);
                duxdy = 2.*x.^2.*(x - 1).^2.*(6.*y.^2 - 6.*y + 1);
                duydx = -2.*y.^2.*(y - 1).^2.*(6.*x.^2 - 6.*x + 1);
                duydy = -duxdx;
                eU = eU + sum(W.*((duxdx - Guh(1,1)).^2 + (duxdy - Guh(2,1)).^2 + ...
                                  (duydx - Guh(1,2)).^2 + (duydy - Guh(2,2)).^2));

                pEx = x.^3 + y.^3 - 0.5;
                eP  = eP + sum(W.*(pEx - L*ph(nodes)).^2);
            end
            errUH1 = sqrt(eU);
            errPL2 = sqrt(eP);
            fprintf('%-32s |u - u_h|_H1 = %10.4e   ||p - p_h||_L2 = %10.4e\n', ...
                    label, errUH1, errPL2);
        end

        function [uPrime, xc] = plotSubscale(obj, type)
            % Subscale velocity u' recovered after the solve, element by element.
            % From the subscale row of the block system (K_vu = 0):
            %   'bubble' / 'exact' :  c = K_vv^-1 (f_v - B_vp p_e)
            %   'min' 'mean' 'max' :  c = tau M_vv^-1 (f_v - B_vp p_e)   (u' = tau x residual)
            % u' is evaluated at the element centroid: u'(xc) = sum_k c_k phi_k(xc).
            % Call it right after MixedFormulation_bubble / _eigen with the same type.
            nu     = obj.nu.constant;
            coord  = obj.mesh.coord;
            connec = obj.mesh.connec;
            nElem  = size(connec, 1);
            ph     = double(obj.pFun.fValues);
            [L, w] = obj.triangleQuadrature(6);
            Lc     = [1 1 1]/3;                              % centroid (barycentric)

            uPrime = zeros(nElem, 2);
            xc     = zeros(nElem, 2);
            for e = 1:nElem
                nodes = connec(e,:);
                xe = coord(nodes,:);
                T  = [ones(3,1), xe];
                C  = T \ eye(3);
                dL = C(2:3,:);
                W  = abs(det(T))*w;
                xg = L*xe;
                [fx, fy] = obj.load(xg(:,1), xg(:,2));

                [phi, dphi] = obj.fineModes(L, dL, type);
                m  = size(phi, 2);
                K1 = zeros(m);
                for d = 1:2
                    K1 = K1 + nu*dphi(:,:,d)'*(W.*dphi(:,:,d));
                end
                M1  = phi'*(W.*phi);
                Bvp = [ (phi'*W)*dL(1,:) ; (phi'*W)*dL(2,:) ];
                fv  = [ phi'*(W.*fx) ; phi'*(W.*fy) ];
                r   = fv - Bvp*ph(nodes);                    % subscale residual

                switch type
                    case {'bubble', 'exact'}
                        c = blkdiag(K1, K1) \ r;
                    otherwise
                        c = obj.tau(e) * (blkdiag(M1, M1) \ r);
                end

                phiC = obj.fineModes(Lc, dL, type);          % modes at the centroid (1 x m)
                uPrime(e,:) = [phiC*c(1:m), phiC*c(m+1:2*m)];
                xc(e,:)     = Lc*xe;
            end

            uMag = sqrt(sum(uPrime.^2, 2));
            fprintf('Subscale (%s): max|u''| = %.4e, mean|u''| = %.4e\n', type, max(uMag), mean(uMag));

            figure('Name', ['Subscale u'' - ' type])
            subplot(1,2,1)
            patch('Faces', connec, 'Vertices', coord, 'FaceVertexCData', uMag, ...
                  'FaceColor', 'flat', 'EdgeColor', [0.6 0.6 0.6]);
            axis equal tight; colorbar; title(['|u''| at the centroids - ' type])
            subplot(1,2,2)
            quiver(xc(:,1), xc(:,2), uPrime(:,1), uPrime(:,2)); hold on
            triplot(connec, coord(:,1), coord(:,2), 'Color', [0.8 0.8 0.8]);
            axis equal tight; title(['u'' vectors - ' type])
        end

        function plotKvvComparison(obj, e, doPlot)
            % Element e, fine space {b l1, b l2, b l3} (one velocity component):
            %   exact K_vv  vs  tau^-1 M_vv  for tau = 1/lambda_min, mean, 1/lambda_max
            % and the pressure block that enters the system:
            %   S_pp = B' K_vv^-1 B  (exact)  vs  tau B' M_vv^-1 B
            % The matrices are always printed; doPlot = true also draws them.
            if nargin < 3; doPlot = false; end
            nu  = obj.nu.constant;
            nodes = obj.mesh.connec(e,:);
            xe  = obj.mesh.coord(nodes,:);
            T   = [ones(3,1), xe];
            C   = T \ eye(3);
            dL  = C(2:3,:);
            A2  = abs(det(T));
            [L, w] = obj.triangleQuadrature(6);
            W   = A2*w;
            [phi, dphi] = obj.fineModes(L, dL, 'exact');

            K1 = zeros(3);
            for d = 1:2
                K1 = K1 + nu*dphi(:,:,d)'*(W.*dphi(:,:,d));
            end
            M1  = phi'*(W.*phi);
            Bvp = [ (phi'*W)*dL(1,:) ; (phi'*W)*dL(2,:) ];  % 6 x 3
            Kb  = blkdiag(K1, K1);
            Mb  = blkdiag(M1, M1);

            lam  = sort(real(eig(K1, M1)));
            taus = [1/lam(1), 2/(lam(1) + lam(end)), 1/lam(end)];
            names = {'1/\lambda_{min}', '2/(\lambda_{min}+\lambda_{max})', '1/\lambda_{max}'};

            Sex = Bvp'*(Kb\Bvp);
            h   = sqrt(A2);                                   % element size (sqrt of 2*area)
            fprintf('\nElement %d: lambda h^2/nu = %s, kappa = %.2f\n', e, ...
                    mat2str(lam'*h^2/nu, 4), lam(end)/lam(1));

            % --- print the matrices ---
            tauNames = {'1/lambda_min', '2/(lambda_min+lambda_max)', '1/lambda_max'};
            fprintf('\nK_vv exact:\n');  disp(K1)
            for i = 1:numel(taus)
                Kt = M1/taus(i);
                fprintf('tau^-1 M_vv, tau = %s   (error %.1f%%):\n', tauNames{i}, ...
                        100*norm(K1 - Kt, 'fro')/norm(K1, 'fro'));
                disp(Kt)
            end
            fprintf('S_pp exact (B'' K_vv^-1 B):\n');  disp(Sex)
            for i = 1:numel(taus)
                St = taus(i)*Bvp'*(Mb\Bvp);
                fprintf('S_pp, tau = %s   (error %.1f%%):\n', tauNames{i}, ...
                        100*norm(Sex - St, 'fro')/norm(Sex, 'fro'));
                disp(St)
            end

            if ~doPlot; return; end

            % --- K_vv ---
            figure('Name', sprintf('K_vv comparison, element %d', e), 'Position', [100 100 1400 650])
            cl = [min(K1(:)), max(K1(:))];
            for i = 1:numel(taus)
                cl = [min(cl(1), min(M1(:))/taus(i)), max(cl(2), max(M1(:))/taus(i))];
            end
            subplot(2,4,1)
            imagesc(K1, cl); axis square; colorbar; title('K_{vv} exact')
            set(gca, 'XTick', 1:3, 'YTick', 1:3)
            for i = 1:numel(taus)
                Kt  = M1/taus(i);
                err = norm(K1 - Kt, 'fro')/norm(K1, 'fro');
                subplot(2,4,1+i)
                imagesc(Kt, cl); axis square; colorbar
                set(gca, 'XTick', 1:3, 'YTick', 1:3)
                title({['\tau^{-1} M_{vv},  \tau = ' names{i}], sprintf('error = %.1f%%', 100*err)})
            end

            % --- S_pp (what enters the system) ---
            sl = [min(Sex(:)), max(Sex(:))];
            Sap = cell(1, numel(taus));
            for i = 1:numel(taus)
                Sap{i} = taus(i)*Bvp'*(Mb\Bvp);
                sl = [min(sl(1), min(Sap{i}(:))), max(sl(2), max(Sap{i}(:)))];
            end
            subplot(2,4,5)
            imagesc(Sex, sl); axis square; colorbar; title('S_{pp} exact (K_{vv}^{-1})')
            set(gca, 'XTick', 1:3, 'YTick', 1:3)
            for i = 1:numel(taus)
                err = norm(Sex - Sap{i}, 'fro')/norm(Sex, 'fro');
                subplot(2,4,5+i)
                imagesc(Sap{i}, sl); axis square; colorbar
                set(gca, 'XTick', 1:3, 'YTick', 1:3)
                title({['S_{pp},  \tau = ' names{i}], sprintf('error = %.1f%%', 100*err)})
                fprintf('  tau = %-24s  error K_vv = %6.1f%%   error S_pp = %6.1f%%\n', ...
                        strrep(strrep(strrep(names{i},'\',''),'_{',''),'}',''), ...
                        100*norm(K1 - M1/taus(i),'fro')/norm(K1,'fro'), 100*err);
            end
        end

        function plotResults(obj, lab)
            X  = obj.mesh.coord;
            cn = obj.mesh.connec;
            uh = double(obj.uFun.fValues);
            ph = double(obj.pFun.fValues);
            pEx = X(:,1).^3 + X(:,2).^3 - 0.5;

            figure('Name', lab)
            subplot(1,3,1)
            trisurf(cn, X(:,1), X(:,2), ph); shading interp; view(2); axis equal tight; colorbar
            title(['Pressure - ' lab])
            subplot(1,3,2)
            trisurf(cn, X(:,1), X(:,2), pEx); shading interp; view(2); axis equal tight; colorbar
            title('Pressure - exact')
            subplot(1,3,3)
            quiver(X(:,1), X(:,2), uh(:,1), uh(:,2)); axis equal tight
            title(['Velocity - ' lab])
        end

    end

end
