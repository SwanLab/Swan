classdef TutorialMixedStokes < handle


    properties (Access = private)
        mesh
        nu
        uFun, pFun
        bcU, bcP
        tau
    end

    methods (Access = public)

        function obj = TutorialMixedStokes()
            obj.createMesh()
            obj.createMaterialProperties()
            obj.createSolutionField()
            obj.createBoundaryConditions()

            obj.MixedFormulation_classic()
            obj.plotResults('u-p (P1-P1) no stab.');

            obj.MixedFormulation_bubble()
            obj.plotResults('u-p (P1-P1) bubble funct.');
            uB = double(obj.uFun.fValues);  pB = double(obj.pFun.fValues);

            obj.MixedFormulation_eigen('min')
            obj.plotResults('u-p (P1-P1) Tau');

            uE = double(obj.uFun.fValues);  pE = double(obj.pFun.fValues);

            fprintf('max|u_bubble - u_tauMin| = %e\n', max(abs(uB - uE)));
            fprintf('max|p_bubble - p_tauMin| = %e\n', max(abs(pB - pE)));

        end

    end

    methods (Access = private)
    
        function createMesh(obj)

            L = 1.0; 
            nx = 50; 
            s.coord = linspace(0, L, nx+1)';
            s.connec = [(1:nx)', (2:nx+1)'];
            obj.mesh = Mesh.create(s);
          
        end

        function createMaterialProperties(obj)
            viscosity = 1e-6;
            obj.nu =ConstantFunction.create(viscosity, obj.mesh);

        end

        function createSolutionField(obj)
            obj.uFun = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.pFun = LagrangianFunction.create(obj.mesh,1,'P1');
        end

        function createBoundaryConditionsVelocity(obj)
            TOL = 1e-12;
            xMax = max(obj.mesh.coord);
            xMin = min(obj.mesh.coord);

            isLeft  = @(coor)  abs(coor(:,1)-xMin)< TOL;
            isRight = @(coor)  abs(coor(:,1)-xMax)< TOL;

            % Boundary conditions for velocity
            sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor);
            sDir{1}.direction = [1];
            sDir{1}.value     = 0;                    
            sDir{1}.ndim = length(sDir{1}.direction);

            obj.bcU = obj.createBoundaryConditionsFromInput(sDir);

        end

        function createBoundaryConditionsPressure(obj)
            TOL = 1e-12;
            xMax = max(obj.mesh.coord);
            xMin = min(obj.mesh.coord);

            isLeft  = @(coor)  abs(coor(:,1)-xMin)< TOL;

            % Boundary conditions for velocity
            sDir{1}.domain    = @(coor) isLeft(coor) ;
            sDir{1}.direction = [1];
            sDir{1}.value     = 0;                    
            sDir{1}.ndim = length(sDir{1}.direction);

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

        function createBoundaryConditions(obj)
            obj.createBoundaryConditionsVelocity();
            obj.createBoundaryConditionsPressure();
        end

        function MixedFormulation_classic(obj)
            LHS = obj.createLHS_v1();
            RHS = obj.createRHS_v1();

            nUt   = obj.uFun.nDofs;
            dofFU = obj.bcU.free_dofs(:);
            dofFP = obj.bcP.free_dofs(:);
            free  = [dofFU; nUt + dofFP];

            LHSr = LHS(free, free);
            RHSr = RHS(free);

            %x = full(LHSr)\RHSr;

            x = pinv(full(LHSr))*RHSr; %Pseudo inverse (with minimal norm)

            nU = length(dofFU);
            nP = length(dofFP);

            uF = x(1:nU);
            pF = x(nU+1:nU+nP);

            uT = zeros(obj.uFun.nDofs,1);
            uT(dofFU,1) = uF; 
            uT = reshape(uT,obj.uFun.ndimf,[])';
            obj.uFun.setFValues(uT);

            pT = zeros(obj.pFun.nDofs,1);
            pT(dofFP,1) = pF;
            obj.pFun.setFValues(pT);

        end

        function MixedFormulation_bubble(obj)

            LHS = obj.createLHS_v1();
            RHS = obj.createRHS_v1();

            [Spp, Fp] = obj.createBubbleCondensation();

            nUt = obj.uFun.nDofs;
            iP  = nUt + (1:obj.pFun.nDofs);
            LHS(iP,iP) = LHS(iP,iP) + Spp;
            RHS(iP)    = RHS(iP)    + Fp;

            nUt   = obj.uFun.nDofs;
            dofFU = obj.bcU.free_dofs(:);
            dofFP = obj.bcP.free_dofs(:);
            free  = [dofFU; nUt + dofFP];

            LHSr = LHS(free, free);
            RHSr = RHS(free);

            x = LHSr\RHSr;

            nU = length(dofFU);
            nP = length(dofFP);

            uT = zeros(obj.uFun.nDofs,1);
            uT(dofFU,1) = x(1:nU);
            obj.uFun.setFValues(uT);

            pT = zeros(obj.pFun.nDofs,1);
            pT(dofFP,1) = x(nU+1:nU+nP);
            obj.pFun.setFValues(pT);
        end

        function MixedFormulation_eigen(obj, tauType)
            % Case 3: fine space = 2 internal modes per element
            %   b2 = xi(1-xi) (quadratic),  b3 = xi(1-xi)(1-2xi) (cubic)
            % K_vv is replaced by tau^-1 M_vv, with tau from the
            % eigenvalues of K_vv x = lambda M_vv x.
            % tauType: 'max' -> tau = 1/lambda_max
            %          'mean'-> tau = 2/(lambda_min + lambda_max)
            %          'min' -> tau = 1/lambda_min
            %          'exact' -> exact condensation (K_vv^-1, no tau)
            LHS = obj.createLHS_v1();
            RHS = obj.createRHS_v1();
 
            [Spp, Fp] = obj.createEigenCondensation(tauType);
 
            nUt = obj.uFun.nDofs;
            iP  = nUt + (1:obj.pFun.nDofs);
            LHS(iP,iP) = LHS(iP,iP) + Spp;
            RHS(iP)    = RHS(iP)    + Fp;
 
            dofFU = obj.bcU.free_dofs(:);
            dofFP = obj.bcP.free_dofs(:);
            free  = [dofFU; nUt + dofFP];
 
            x = LHS(free,free)\RHS(free);
 
            nU = length(dofFU);
            nP = length(dofFP);
 
            uT = zeros(obj.uFun.nDofs,1);
            uT(dofFU,1) = x(1:nU);
            obj.uFun.setFValues(uT);
 
            pT = zeros(obj.pFun.nDofs,1);
            pT(dofFP,1) = x(nU+1:nU+nP);
            obj.pFun.setFValues(pT);
        end

        function LHS = createLHS_v1(obj)
    
            nu = obj.nu;
    
            K = IntegrateLHS(@(u,v) nu.*DP(Grad(v),Grad(u)), obj.uFun, obj.uFun, obj.mesh, 'Domain', 2);
    
            B = IntegrateLHS(@(p,v) (-1).*DP(Grad(v),p), obj.uFun, obj.pFun, obj.mesh, 'Domain', 2);
    
            Z = sparse(obj.pFun.nDofs, obj.pFun.nDofs);
    
            LHS = [K, B; B', Z];
    
        end
    
        function RHS = createRHS_v1(obj)
    
            f  = AnalyticalFunction.create(@(x) pi*cos(pi*x(1,:,:)), obj.mesh);
            Fu = IntegrateRHS(@(v) DP(v,f), obj.uFun, obj.mesh, 'Domain', 4);
            Fp = zeros(obj.pFun.nDofs, 1);
    
            RHS = [Fu; Fp];
        end

        function [Spp, Fp] = createBubbleCondensation(obj)
            % For each element, bubble b = 4 xi (1 - xi), xi in [0,1]:
            %   Kvv = nu (b', b')            (1 x 1)
            %   Bvp = -(N_j, b')             (1 x 2)  row v, pressure columns
            %   fv  = (f, b)                 (1 x 1)
            %   Kuv = nu (b', N_i') = 0 in 1D
            % Condensation:  Spp = -Bvp' Kvv^-1 Bvp,  Fp = -Bvp' Kvv^-1 fv
            nu    = obj.nu.constant;
            f     = @(x) pi*cos(pi*x);
            coord = obj.mesh.coord(:,1);
            connec = obj.mesh.connec;
            nElem = size(connec,1);
            nP    = obj.pFun.nDofs;

            t5 = [-0.9061798459386640; -0.5384693101056831; 0; 0.5384693101056831; 0.9061798459386640];
            w5 = [0.2369268850561891; 0.4786286704993665; 0.5688888888888889; 0.4786286704993665; 0.2369268850561891];
            xi = (t5 + 1)/2;                                    % 5 Gauss points in [0,1]
            w  = w5/2;                                 % Gauss weights
            N  = [1 - xi, xi];                                  % P1 shape functions
            b  = 4*xi.*(1 - xi);                                % bubble

            Spp = sparse(nP, nP);
            Fp  = zeros(nP, 1);
            for e = 1:nElem
                nodes = connec(e,:);
                xe = coord(nodes);
                h  = xe(2) - xe(1);
                db = 4*(1 - 2*xi)/h;                            % b'
                xg = N*xe;                                      % Gauss points in x

                Kvv = h*sum(w.*(nu*db.^2));
                Bvp = -h*(w.*db)'*N;
                fv  = h*sum(w.*(f(xg).*b));

                Spp(nodes,nodes) = Spp(nodes,nodes) - Bvp'*Bvp/Kvv;
                Fp(nodes)        = Fp(nodes)        - Bvp'*fv/Kvv;
            end
        end

        function [Spp, Fp] = createEigenCondensation(obj, tauType)
            nu     = obj.nu.constant;              % use obj.nu.constant if nu is a ConstantFunction
            f      = @(x) pi*cos(pi*x);
            coord  = obj.mesh.coord(:,1);
            connec = obj.mesh.connec;
            nElem  = size(connec,1);
            nP     = obj.pFun.nDofs;
 
            t5 = [-0.9061798459386640; -0.5384693101056831; 0; 0.5384693101056831; 0.9061798459386640];
            w5 = [0.2369268850561891; 0.4786286704993665; 0.5688888888888889; 0.4786286704993665; 0.2369268850561891];
            xi = (t5 + 1)/2;                                    % 5 Gauss points in [0,1]
            w  = w5/2;                                          % (exact up to degree 9)
            N  = [1 - xi, xi];                                  % P1 shape functions
            b  = [xi.*(1-xi), xi.*(1-xi).*(1-2*xi)];            % fine modes (3 x 2)
            dbRef = [1 - 2*xi, 1 - 6*xi + 6*xi.^2];             % d b / d xi
 
            Spp = sparse(nP, nP);
            Fp  = zeros(nP, 1);
            obj.tau = zeros(nElem, 1);
            for e = 1:nElem
                nodes = connec(e,:);
                xe = coord(nodes);
                h  = xe(2) - xe(1);
                db = dbRef/h;                                   % b'
                xg = N*xe;
 
                Kvv = h*db'*(w.*db)*nu;                         % nu (b_i', b_j')  2x2
                Mvv = h*b'*(w.*b);                              % (b_i, b_j)       2x2
                Bvp = -h*db'*(w.*N);                            % -(N_j, b_i')     2x2
                fv  = h*b'*(w.*f(xg));                          % (f, b_i)         2x1
 
                lam = eig(Kvv, Mvv);
                switch tauType
                    case 'max';  t = 1/max(lam);
                    case 'mean'; t = 2/(min(lam) + max(lam));
                    case 'min';  t = 1/min(lam);
                end
 
                if strcmp(tauType, 'exact')
                    Se = -Bvp'*(Kvv\Bvp);
                    Fe = -Bvp'*(Kvv\fv);
                else
                    Se = -t*Bvp'*(Mvv\Bvp);                     % K_vv ~ tau^-1 M_vv
                    Fe = -t*Bvp'*(Mvv\fv);
                    obj.tau(e) = t;
                end
 
                Spp(nodes,nodes) = Spp(nodes,nodes) + Se;
                Fp(nodes)        = Fp(nodes)        + Fe;
            end
        end

    
        function plotResults(obj, lab)

            x = obj.mesh.coord(:,1);
        
            figure('Name', lab)
            subplot(2,1,1)
            plot(x, zeros(size(x)), 'k-', 'LineWidth', 1.5); hold on
            plot(x, obj.uFun.fValues, 'o-'); grid on
            xlabel('x'); ylabel('u'); title(['Velocity'])
            legend('Exact', lab, 'Location', 'best')
        
            subplot(2,1,2)
            plot(x, sin(pi*x), 'k-', 'LineWidth', 1.5); hold on
            plot(x, obj.pFun.fValues, 'o-'); grid on
            xlabel('x'); ylabel('p'); title(['Pressure'])
            legend('Exact', lab, 'Location', 'south')
        end


    end

end