classdef InteriorPointMethodsSolver < handle 
    properties (Access = public)
        bp
        invDiagdL, invDiagdU
        sigmaL, sigmaU
        H
        alpha_du, alpha_du_max
        alpha_pr, alpha_pr_max
        firstLHS
        iteration
        sol
        x, s
        xL, sL, bL
        xU, sU, bU
        gradient, jacobian, residualLS
        lambda
        thetaMax
        filter, store
        e
        th
        ph
        status
        dx, ds 
        dlam, lam
        dzL, dzU
        zLa, zUa
        xa, sa
        m, n, ns
        tau
        zL, zU
    end
    properties (Access = protected)
        okJacobian
        okHessian
        okObjJacobian
        thetaMin
        initResidual
    end
    properties (Access = private)
        isposdef
        e_mu
        diagonalzL
        diagonalzU
        diagonaldL
        diagonaldU
        hessian
        dL
        dU
        search1
        search2
        linearLHS
        updatedLambda
    end

    methods (Access = public)
        function obj = InteriorPointMethodsSolver(cParams)
            obj.init(cParams);
        end
        function solve(obj)
            obj.openOptimizer();
        end

        function iterationProcess(obj)
            for iter = 1:obj.bp.maxiter
                obj.iteration = iter;
                u.bp = obj.bp;
                u.x = obj.x;
                u.s = obj.s;
                u.bL = obj.bL;
                u.bU = obj.bU;
                residual = obj.residualComputer(u);
                obj.residualLS = residual;
                obj.th = obj.thetaComputer(u);
                u.xL = obj.xL;
                u.xU = obj.xU;
                obj.ph = obj.phiComputer(u);
                obj.zUpdater();
                obj.sigmaComputer();
                obj.computeJacobian();
                obj.computeHessian();
                obj.computeObjectiveGradient();
                obj.H = obj.hessian + obj.sigmaL + obj.sigmaU;
                obj.linearSystemSolver(residual);
                obj.computeAcceptancePoint();
                obj.checkConstraintViolations();
                obj.lineSearch();
                if (obj.e_mu <= obj.bp.e_tol)
                    fprintf(1,'\nSuccessful solution\n');
                    obj.status = 'success';
                    obj.finalSolTest();
                    break;
                end
            end
            obj.displaySolution();
        end
    end

    methods(Access = private)
        function zUpdater(obj)
            if (obj.bp.z_update == 2)
                zLC = obj.zL;
                zUC = obj.zU;
                % update explicitly from z = mu / x
                for i = 1:obj.n
                    zLC(i) = obj.bp.mu / (obj.x(i) - obj.xL(i));
                end
                for i = 1:obj.ns
                    zLC(obj.n+i) = obj.bp.mu / (obj.s(i) - obj.sL(i));
                end
                % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
                for i = 1:obj.n
                    zUC(i) = obj.bp.mu / (obj.xU(i) - obj.x(i));
                end
                for i = 1:obj.ns
                    zUC(obj.n+i) = obj.bp.mu / (obj.sU(i) - obj.s(i));
                end
                obj.zL = zLC;
                obj.zU = zUC;
            end
            % make diagonal matrices
            obj.diagonalzL = diag(obj.zL);
            obj.diagonalzU = diag(obj.zU);
        end

        function sigmaComputer(obj)
            % sigmas
            obj.dL = [obj.x-obj.xL obj.s-obj.sL];
            obj.dU = [obj.xU-obj.x obj.sU-obj.s];
            obj.diagonaldL = diag(obj.dL);
            obj.diagonaldU = diag(obj.dU);
            invdML = zeros(obj.n+obj.ns,obj.n+obj.ns);
            invdMU = zeros(obj.n+obj.ns,obj.n+obj.ns);
            sigL = zeros(obj.n+obj.ns,obj.n+obj.ns);
            sigU = zeros(obj.n+obj.ns,obj.n+obj.ns);
            for i = 1:obj.n+obj.ns
                invdML(i,i) = 1 / obj.diagonaldL(i,i); 
                invdMU(i,i) = 1 / obj.diagonaldU(i,i); 
                sigL(i,i) = obj.diagonalzL(i,i) / obj.diagonaldL(i,i);
                sigU(i,i) = obj.diagonalzU(i,i) / obj.diagonaldU(i,i);
             end
             obj.invDiagdL = invdML;
             obj.invDiagdU = invdMU;
             obj.sigmaL = sigL;
             obj.sigmaU = sigU;
        end

        function linearSystemSolver(obj,residual)
            % Construct and solve linear system Ax=b    
            % Symmetric, zL/zU explicitly solved
            % construct A and b
            A1 = [obj.H,obj.jacobian';obj.jacobian,zeros(obj.m,obj.m)];
            %b1(1:n+ns,1) = g' - zL' + zU' + J'*lam';
            b1(1:obj.n + obj.ns,1) = obj.gradient' - obj.zL' + obj.zU' + obj.jacobian'*obj.lambda' + obj.zL'-obj.bp.mu*obj.invDiagdL*obj.e -obj.zU'+obj.bp.mu*obj.invDiagdU*obj.e;
            b1(obj.n + obj.ns + 1:obj.n + obj.ns + obj.m,1) = residual(1:obj.m)';
            % compute search direction, solving Ax = b
            d1 = -pinv(full(A1))*b1;
            
            dx1 = d1(1:obj.n,1);
            ds1 = d1(obj.n+1:obj.n+obj.ns,1);
            dlam1 = d1(obj.n+obj.ns+1:obj.n+obj.ns+obj.m,1);
            % compute search direction for z (explicit solution)
            dzL1 = obj.bp.mu*obj.invDiagdL*obj.e - obj.zL' - obj.sigmaL*[dx1; ds1];
            dzU1 = obj.bp.mu*obj.invDiagdU*obj.e - obj.zU' + obj.sigmaU*[dx1; ds1];
            obj.search1 = [d1;dzL1;dzU1];
            obj.linearLHS = A1;
            obj.testLinearSystemDefinitiness();
            % Un-symmetric, zL/zU implicitly solved
            % construct A
            A2 = [obj.hessian obj.jacobian' -eye(obj.n+obj.ns) eye(obj.n+obj.ns);
                    obj.jacobian zeros(obj.m,obj.m+2*(obj.n+obj.ns));
                    diag(obj.zL) zeros((obj.n+obj.ns),obj.m) obj.diagonaldL zeros(obj.n+obj.ns);
                    -diag(obj.zU) zeros(obj.n+obj.ns,obj.m+obj.n+obj.ns) obj.diagonaldU];
            b2(1:obj.n+obj.ns,1) = obj.gradient' - obj.zL' + obj.zU' + obj.jacobian'*obj.lambda';
            b2(obj.n+obj.ns+1:obj.n+obj.ns+obj.m,1) = residual(1:obj.m)';
            b2(obj.n+obj.ns+obj.m+1:2*(obj.n+obj.ns)+obj.m) = obj.diagonaldL*diag(obj.zL)*obj.e - obj.bp.mu*obj.e;
            b2(2*(obj.n+obj.ns)+obj.m+1:3*(obj.n+obj.ns)+obj.m) = obj.diagonaldU*diag(obj.zU)*obj.e - obj.bp.mu*obj.e;
            d2 = -pinv(full(A2))*b2;
            dx2 = d2(1:obj.n,1);
            ds2 = d2(obj.n+1:obj.n+obj.ns,1);
            dlam2 = d2(obj.n+obj.ns+1:obj.n+obj.ns+obj.m,1);
            dzL2 = d2(obj.n+obj.ns+obj.m+1:2*(obj.n+obj.ns)+obj.m);
            dzU2 = d2(2*(obj.n+obj.ns)+obj.m+1:3*(obj.n+obj.ns)+obj.m);
            obj.search2 = d2;
            obj.isposdef = true;
            obj.linearLHS = A2;
            idebug = obj.bp.idebug;
            matrix = obj.bp.matrix;
            obj.testLinearSystemDefinitiness();
            obj.conditionPrinter(idebug,A1,A2,obj.search1,obj.search2);
            obj.inverseMatrices(matrix,dx1,ds1,dlam1,dzL1,dzU1,idebug,A1,dx2,ds2,dlam2,dzL2,dzU2,A2);
            obj.firstLHS = A1;
        end

        function testLinearSystemDefinitiness(obj)
                if (obj.bp.idebug>=2)
                    if (min(eigs(obj.linearLHS))<0)
                        fprintf(1,'%12.4e not pos definite - eigenvalue test\n',obj.linearLHS);
                    end
                end
                for i=1:length(obj.linearLHS)
                    if ( det( obj.linearLHS(1:i, 1:i) ) <= 0 )
                        obj.isposdef = false;
                        break;
                    end
                end
                if(~obj.isposdef && obj.bp.idebug>=2)
                    fprintf(1,'%12.4e not pos definite - determinant test\n',obj.linearLHS);
                end    
        end

        function computeAcceptancePoint(obj)
            xaC = obj.x + obj.alpha_pr * obj.dx';
            if (obj.ns>=1)
                saC = obj.s + obj.alpha_pr * obj.ds';
            else
               saC = [];
            end
            obj.lam = obj.lambda;
            obj.lambda = obj.lambda + obj.alpha_du * obj.dlam';
            switch(obj.bp.z_update)
                case(1)
                    obj.zLa = obj.zL + obj.alpha_du * obj.dzL';
                    obj.zUa = obj.zU + obj.alpha_du * obj.dzU';
                case(2)
                    for i = 1:obj.n
                        zLaC(i) = obj.bp.mu / (xaC(i) - obj.xL(i));
                        dzLC(i,1) = zLaC(i) - obj.zL(i);
                    end
                    for i = 1:obj.ns
                        zLaC(obj.n+i) = obj.bp.mu / (saC(i) - obj.sL(i));
                        dzLC(obj.n+i,1) = zLaC(obj.n+i) - obj.zL(obj.n+i);
                    end
                    for i = 1:obj.n
                        zUaC(i) = obj.bp.mu / (obj.xU(i)-xaC(i));
                        dzUC(i,1) = zUaC(i) - obj.zU(i);
                    end
                    for i = 1:obj.ns
                        zUaC(obj.n+i) = obj.bp.mu / (obj.sU(i)-saC(i));
                        dzUC(obj.n+i,1) = zUaC(obj.n+i) - obj.zU(obj.n+i);
                    end
                obj.zLa = zLaC;
                obj.dzL = dzLC;
                obj.zUa = zUaC;
                obj.dzU = dzUC;
            end
                obj.xa = xaC;
                obj.sa = saC;
                obj.alpha_pr_max = 1.0;
                obj.alpha_du_max = 1.0;
        end

        function checkConstraintViolations(obj)
            for i = 1:obj.n
                if(obj.xa(i) < obj.xL(i))
                    obj.alpha_pr_max = min(obj.alpha_pr_max,(obj.xL(i) + obj.tau*(obj.x(i) - obj.xL(i)) - obj.x(i))/obj.dx(i,1));
                end
                if(obj.xa(i) > obj.xU(i))
                    obj.alpha_pr_max = min(obj.alpha_pr_max,(obj.xU(i) - obj.tau*(obj.xU(i) - obj.x(i)) - obj.x(i))/obj.dx(i,1));
                end
            end
            for i = 1:obj.ns
                if(obj.sa(i) < obj.sL(i))
                    obj.alpha_pr_max = min(obj.alpha_pr_max,(obj.sL(i) + obj.tau*(obj.s(i) - obj.sL(i)) - obj.s(i))/obj.ds(i,1));
                end
                if(obj.sa(i) > obj.sU(i))
                    obj.alpha_pr_max = min(obj.alpha_pr_max,(obj.sU(i) + obj.tau*(obj.sU(i) - obj.s(i)) - obj.s(i))/obj.ds(i,1));
                end
            end
            for i = 1:obj.n
                if (obj.bp.z_update==1)          
                    if(obj.zLa(i) < 0)
                        obj.alpha_du_max = min(obj.alpha_du_max,(obj.tau*obj.zL(i) - obj.zL(i))/obj.dzL(i,1));
                    end
                    if(obj.zUa(i) < 0)
                        obj.alpha_du_max = min(obj.alpha_du_max,(obj.tau*obj.zU(i) - obj.zU(i))/obj.dzU(i,1));
                    end
                end
            end
        end

        function lineSearch(obj)
            LSearch = LineSearchComputer(obj);
            LSearch.search();
            obj.x = LSearch.x;
            obj.status = LSearch.status;
            obj.s = LSearch.s;
            obj.lambda = LSearch.lambda;
            obj.zL = LSearch.zL;
            obj.zU = LSearch.zU;
            obj.store = LSearch.store;
            obj.bp.mu = LSearch.bp.mu;
            obj.e_mu = LSearch.e_mu;
        end
        function inverseMatrices(obj,matrix,dx1,ds1,dlam1,dzL1,dzU1,idebug,A1,dx2,ds2,dlam2,dzL2,dzU2,A2)
            if (matrix == 1)
                obj.dx = dx1;
                obj.ds = ds1;
                obj.dlam = dlam1;
                obj.dzL = dzL1;
                obj.dzU = dzU1;
                if(cond(full(A1))>1e10 && idebug>=2)
                    fprintf(1,'Warning: A1 condition number high %12.4e\n',cond(full(A1)))
                end
            else
                obj.dx = dx2;
                obj.ds = ds2;
                obj.dlam = dlam2;
                obj.dzL = dzL2;
                obj.dzU = dzU2;
                if(cond(full(A2))>1e10 && idebug>=2)
                    fprintf(1,'Warning: A2 condition number high %12.4e\n',cond(full(A2)))
                end
            end
        end
        function init(obj,cParams)
            obj.bp = cParams;
        end

        function openOptimizer(obj)
                u = obj.bp;
                opt = OptimizerComputer(u);
                opt.compute();
        end

        function finalSolTest(obj)
            load("final_sol.mat",'x_old');
            u.loadedData = x_old;
            u.actualData = obj.x;
            u.desiredTest = 'Final Solution';
            obj.testResults(u);
        end

        function displaySolution(obj)
            %% Display final solution
            fprintf('Status: '); fprintf(obj.status); fprintf('\n')
            disp('Solution: ')
            disp(obj.x)
            disp('Slack Variables: ')
            disp(obj.s)
            disp('Equation Multipliers: ')
            disp(obj.lambda)
            disp('Lower Constraint Variable Mult: ')
            disp(obj.zL)
            disp('Upper Constraint Variable Mult: ')
            disp(obj.zU)

            if (obj.bp.prob >= 1)
                %% Display figures on iteration progress
                figure(1)
                hold off;
                i = 1;
                plot(obj.store(:,1),obj.store(:,i+1),'k-')
                hold on;
                i = 2;
                plot(obj.store(:,1),obj.store(:,i+1),'b-')
                legend('x_1','x_2')
                
                figure(2)
                hold off;
                plot(obj.store(:,1),obj.store(:,obj.n+2),'r-')
                hold on;
                plot(obj.store(:,1),obj.store(:,obj.n+3),'b-')
                plot(obj.store(:,1),log10(obj.store(:,obj.n+4)),'g-')
                legend('alpha_{pr}','alpha_{du}','log_{10}(mu)');
            end
            %record solution for function return
            obj.sol.status = obj.status;
            obj.sol.x = obj.x;
            obj.sol.s = obj.s;
            obj.sol.lam = obj.lambda;
            obj.sol.zL = obj.zL;
            obj.sol.zU = obj.zU;
        end
    end

    methods (Static, Access = public)
        function residual = residualComputer(cParams)
            residualC = ResidualComputer(cParams);
            residualC.compute();
            residual = residualC.c;
        end
        
        function phi = phiComputer(cParams)
            phiC = PhiComputer(cParams);
            phiC.compute();
            phi = phiC.phi;
        end
    end
    methods (Static, Access = protected)
        function conditionPrinter(idebug,A1,A2,search1,search2)
            if (idebug>=2)
                fprintf(1,'Diff: %12.4e, Cond(A1): %12.4e, Cond(A2): %12.4e\n', ...
                    sum(abs(search1 - search2)),cond(full(A1)),cond(full(A2)));
            end
        end

        function theta = thetaComputer(cParams)
            th = ThetaComputer(cParams);
            th.compute();
            theta = th.theta;
        end

        function testResults(cParams)
            Test = TestComputer(cParams);
            Test.compute();
        end
    end
    methods (Access = protected)
        function computeObjectiveGradient(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s =obj.s;
            grad = GradientComputer(u);
            grad.create();
            obj.gradient = grad.objGradient;
        end

        function computeJacobian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.bL = obj.bL;
            u.bU = obj.bU;
            jac = JacobianComputer(u);
            jac.compute();
            obj.jacobian = jac.pd;
        end

        function computeHessian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.lam = obj.lambda;
            hes = HessianComputer(u);
            hes.create();
            obj.hessian = hes.hess;
        end
    end
end