classdef LineSearchComputer < handle

    properties (Access = public)
        e_mu, e
        x, s
        status
        lambda
        zL, zU
        store
        bp        
    end
    properties (Access = private)
        predictReduction, predictUpdated
        reduced
        alpha_du, alpha_du_max
        alpha_pr, alpha_pr_max
        xa, sa
        dx
        n, ns, m
        ds
        dzL, zLa
        dzU, zUa
        dlam, lam
        xL, sL
        tau
        xU, sU
        bU, bL
        accept
        gradient, jacobian, residualLS
        H
        filter
        invDiagdU, invDiagdL
        th, ph, thetaMax
        sigmaL, sigmaU
        iter
        firstLHS
    end
    methods (Access = public)
        function obj = LineSearchComputer(cParams)
            obj.init(cParams);
        end
        function search(obj)
            obj.createLineSearchMethod();
            obj.updateVariables();
            obj.updateZvalues();
            obj.pushVariables();
            obj.reduceMeritFunction();
            obj.performLineSearch();
            obj.contourPlot();
            obj.checkAcceptance();
            obj.updateWithAcceptance();
            obj.checkConvergence();
            obj.checkNewBarrierProblem();
            obj.printAndStoreVariables();
            obj.checkMaxIter();
            obj.checkNotAcceptance();
        end
    end

    methods (Access = private)
        function createLineSearchMethod(obj)
            switch obj.bp.line_search
                %  1 = reduction in merit function
                %  2 = simple clipping
                %  3 = filter method
                case(1)
                    % alpha_du set as approach to constraint
                    obj.alpha_du = obj.alpha_du_max;
                    % compute new acceptance point
                    obj.alpha_pr = min(obj.alpha_pr,obj.alpha_pr_max);
                case(2)
                    % test case for alpha_pr, alpha_du = 1 with clipping only
                    obj.alpha_pr = 1;
                    obj.alpha_du = 1;
                case(3)
                    % alpha_du set as approach to constraint
                    obj.alpha_du = obj.alpha_du_max;
                    % compute new acceptance point
                    obj.alpha_pr = min(obj.alpha_pr,obj.alpha_pr_max);
             end
        end
        
        function updateVariables(obj)
            % apply alpha values
            obj.xa = obj.x + obj.alpha_pr * obj.dx';
                if(obj.ns>=1)
                    obj.sa = obj.s + obj.alpha_pr * obj.ds';
                end
            obj.lambda = obj.lam + obj.alpha_pr * obj.dlam';
        end

        function updateZvalues(obj)
            % updating zLa and zUa
            switch(obj.bp.z_update)
            case(1)
                % update from direct solve approach
                obj.zLa = obj.zL + obj.alpha_du * obj.dzL';
                obj.zUa = obj.zU + obj.alpha_du * obj.dzU';
            case(2)
                % update explicitly from z = mu / x
                for i = 1:obj.n
                    obj.zLa(i) = obj.bp.mu / (obj.xa(i) - obj.xL(i));
                    obj.dzL(i,1) = obj.zLa(i) - obj.zL(i);
                end
                for i = 1:obj.ns
                    obj.zLa(obj.n + i) = obj.bp.mu / (obj.sa(i) - obj.sL(i));
                    obj.dzL(obj.n + i,1) = obj.zLa(obj.n + i) - obj.zL(obj.n + i);
                end
                % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
                for i = 1:obj.n
                    obj.zUa(i) = obj.bp.mu / (obj.xU(i) - obj.xa(i));
                    obj.dzU(i,1) = obj.zUa(i) - obj.zU(i);
                end
                for i = 1:obj.ns
                    obj.zUa(obj.n + i) = obj.bp.mu / (obj.sU(i) - obj.sa(i));
                    obj.dzU(obj.n + i,1) = obj.zUa(obj.n + i) - obj.zU(obj.n + i);
                end
            end
        end

        function pushVariables(obj)
            % clipping (this should already be arranged by alpha_max)
            % push away from the boundary with tau
            for i = 1:obj.n
                if(obj.xa(i) < obj.xL(i))
                    obj.xa(i) = obj.xL(i) + obj.tau*(obj.x(i) - obj.xL(i));
                end
                if(obj.xa(i) > obj.xU(i))
                    obj.xa(i) = obj.xU(i) - obj.tau*(obj.xU(i) - obj.x(i));
                end
                if(obj.zLa(i) < 0)
                    obj.zLa(i) = obj.tau*(obj.zL(i));
                end
                if(obj.zUa(i) < 0)
                    obj.zUa(i) = obj.tau*(obj.zU(i));
                end
            end
            for i = 1:obj.ns
                if(obj.sa(i) < obj.sL(i))
                    obj.sa(i) = obj.sL(i) + obj.tau*(obj.s(i) - obj.sL(i));
                end
                if(obj.sa(i) > obj.sU(i))
                    obj.sa(i) = obj.sU(i) - obj.tau*(obj.sU(i) - obj.s(i));
                end
            end
        end

        function reduceMeritFunction(obj)
            % predicted reduction in the merit function
            obj.predictReduction = -obj.alpha_pr*obj.gradient*[obj.dx;obj.ds] - 0.5*obj.alpha_pr^2*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds] + obj.bp.nu*(norm(obj.residualLS',1) - norm(obj.residualLS' + obj.alpha_pr*obj.jacobian*[obj.dx;obj.ds],1));
            u.bp = obj.bp;
            u.x = obj.x;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.s = obj.s;
            u.bU = obj.bU;
            u.bL = obj.bL;
            meritX = obj.meritFunctionComputer(u);
            u.x = obj.xa;
            u.s = obj.sa;
            meritA = obj.meritFunctionComputer(u);
            obj.reduced = meritX - meritA;
        end

        function performLineSearch(obj)
            % line search criteria
            %  1 = reduction in merit function
            %  2 = simple clipping
            %  3 = filter method
            switch obj.bp.line_search
            case(1)
                pred_phi_1 = obj.gradient*[obj.dx;obj.ds];
                pred_phi_2 = max(0,0.5*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds]);
                pred_phi_decrease = pred_phi_1 + pred_phi_2;
                u.bp = obj.bp;
                u.xa = obj.xa;
                u.sa = obj.sa;
                u.bL = obj.bL;
                u.bU = obj.bU;
                theta = thetaComputer(u);
                rho = 0.1;
                nuUpdated = pred_phi_decrease / ((1 - rho) * theta);
                obj.bp.nu = max(1,min(1000,nuUpdated));
                obj.predictUpdated = -obj.alpha_pr*obj.gradient*[obj.dx;obj.ds] - 0.5*obj.alpha_pr^2*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds] + obj.bp.nu*(norm(obj.residualLS',1) - norm(obj.residualLS' + obj.alpha_pr*obj.jacobian*[obj.dx;obj.ds],1));
                obj.reduceMeritFunction()
                eta = 0.2;
                if (obj.reduced >= eta*obj.predictUpdated)
                    obj.accept = true;
                else
                    obj.accept = false;
                end
            case(2)
                obj.accept = true;
                obj.alpha_pr = 1;
                obj.alpha_du = 1;
            case(3)
                u.bp = obj.bp;
                u.x = obj.x;
                u.xa = obj.xa;
                u.xL = obj.xL;
                u.xU = obj.xU;
                u.filter = obj.filter;
                
                ac = bp_accept(u);
                ac.compute();
                obj.accept = ac.accept;
            end
        end

        function contourPlot(obj)
            if (obj.bp.contour)
                bp_contour(obj.bp,obj.x,obj.alpha_pr*obj.dx);
            end
        end

        function checkAcceptance(obj)
            if (false)
                if (not(obj.accept))
                    u.bp = obj.bp;
                    u.x = obj.x;
                    u.s = obj.s;
                    u.bL = obj.bL;
                    u.bU = obj.bU;
                    residaulX = residualComputer(u);
                    u.x = obj.xa;
                    u.s = obj.sa;
                    residualA = residualComputer(u);
                    b1(obj.n+1:obj.n+obj.m,1) = obj.alpha_pr*residualX' + residualA';
                    d3 = -pinv(full(obj.firstLHS))*b1;
                    dx3 = d3(1:obj.n,1);
                    ds3 = d3(obj.n+1:obj.n+obj.ns,1);
                    dlam3 = d3(obj.n + obj.ns + 1:obj.n + obj.ns + obj.m,1);
                    dzL3 = obj.bp.mu*obj.invDiagdL*obj.e - obj.zL' - obj.sigmaL*[dx3; ds3];
                    dzU3 = obj.bp.mu*obj.invDiagdU*obj.e - obj.zU' + obj.sigmaU*[dx3; ds3];
                    obj.xa = obj.x + obj.alpha_pr * dx3';
                    obj.sa = obj.s + obj.alpha_pr * ds3';
                    obj.lambda = obj.lambda + obj.alpha_du * dlam3';
                    obj.zLa = obj.zL + obj.alpha_du * dzL3';
                    obj.zUa = obj.zU + obj.alpha_du * dzU3';
                    u.bp = obj.bp;
                    u.x = obj.x;
                    u.xa = obj.xa;
                    u.xL = obj.xL;
                    u.xU = obj.xU;
                    u.filter = obj.filter;
                    ac = AcceptanceComputer(u);
                    ac.compute();
                    obj.accept = ac.accept;
                    if (obj.accept)
                        fprintf(1,'2nd order correction success\n');
                    else
                        fprintf(1,'2nd order correction failed to find accepatable point\n');                        
                    end
                end
            end
        end

        function updateWithAcceptance(obj)
            if (obj.accept)
                obj.x = obj.xa;
                obj.s = obj.sa;
                obj.lam = obj.lambda;
                obj.zL = obj.zLa;
                obj.zU = obj.zUa;
                obj.filter = [obj.filter; obj.th obj.ph];
            end
        end

        function checkConvergence(obj)
                s_max = 100; % > 1
                s_d = max(s_max,(sum(abs(obj.lambda))+sum(abs(obj.zL))+sum(abs(obj.zU)))/(obj.m+2*(obj.n+obj.ns)));
                s_c = max(s_max,(sum(abs(obj.zU))+sum(abs(obj.zL)))/(2*(obj.n+obj.ns)));
                obj.computeObjectiveGradient();
                obj.computeJacobian();
                part(1) = max(abs(obj.gradient' + obj.jacobian'*obj.lambda' - obj.zL' + obj.zU'))/s_d;
                u.bp = obj.bp;
                u.x = obj.x;
                u.s = obj.s;
                u.bL = obj.bL;
                u.bU = obj.bU;
                residualUpdated = InteriorPointMethodsSolver.residualComputer(u);
                part(2) = max(abs(residualUpdated));
                part(3) = max(abs(diag([obj.x-obj.xL obj.s-obj.sL])*diag(obj.zL)*obj.e))/s_c;
                part(4) = max(abs(diag([obj.xU-obj.x obj.sU-obj.s])*diag(obj.zU)*obj.e))/s_c;
                obj.e_mu = max(part);
                if (obj.bp.idebug>=1)
                    fprintf(1,'pred: %12.4e, ared: %12.4e, e_mu: %12.4e, k*mu: %12.4e\n',obj.predictUpdated,obj.reduced,obj.e_mu,0.2*obj.bp.mu);
                end
        end

        function checkNewBarrierProblem(obj)
            k_mu  = 0.2;
            if (obj.bp.mu_update)
                if (obj.e_mu < k_mu * obj.bp.mu)
                    th_mu = 1.5;
                    obj.bp.mu = max(obj.bp.e_tol/10,min(k_mu*obj.bp.mu,obj.bp.mu^th_mu));
                    obj.tau = min(obj.bp.tau_max,100*obj.bp.mu);
                    u.bp = obj.bp;
                    u.x = obj.x;
                    u.xL = obj.xL;
                    u.xU = obj.xU;
                    u.s = obj.s;
                    u.bL = obj.bL;
                    u.bU = obj.bU;
                    phi = InteriorPointMethodsSolver.phiComputer(u);
                    obj.filter = [obj.thetaMax phi];
                end
            end
        end

        function printAndStoreVariables(obj)
            obj.printIteration();   
            obj.store = [obj.store; obj.iter obj.x obj.alpha_pr obj.alpha_du obj.bp.mu];
        end

        function printIteration(obj)
            u.bp = obj.bp;
            u.iter = obj.iter;
            u.x = obj.x;
            u.lam = obj.lambda;
            u.zL = obj.zL;
            u.zU = obj.zU;
            u.alpha_pr = obj.alpha_pr;
            u.alpha_du = obj.alpha_du;
            u.s = obj.s;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.bL = obj.bL;
            u.bU = obj.bU;
            pr = PrinterComputer(u);
            pr.print();
        end

        function checkMaxIter(obj)
            if (obj.iter == obj.bp.maxiter)
                obj.status = 'failed: max iterations';
                return;
            end 
        end
        
        function checkNotAcceptance(obj)
            if (obj.accept)
                obj.alpha_pr = 1.0;
                obj.alpha_du = 1.0;
            else
                obj.alpha_pr = obj.alpha_pr / 2;
                obj.alpha_du = obj.alpha_du / 2;
                if (obj.alpha_pr < 1e-4)
                    obj.alpha_pr = 1.0;
                    obj.alpha_du = 1.0;
                end
            end
        end

        function init(obj,cParams)
            obj.bp = cParams.bp;      
            obj.alpha_du = cParams.alpha_du; obj.alpha_du_max = obj.alpha_du_max;
            obj.alpha_pr = cParams.alpha_pr; obj.alpha_pr_max = obj.alpha_pr_max;
            obj.xa = cParams.xa; obj.sa = cParams.sa;
            obj.x = cParams.x; obj.dx = cParams.dx;
            obj.n = cParams.n; obj.ns = cParams.ns; obj.m = cParams.m;
            obj.s = cParams.s; obj.ds = cParams.ds;
            obj.zL = cParams.zL; obj.dzL = cParams.dzL; obj.zLa = cParams.zLa;
            obj.zU = cParams.zU; obj.dzU = cParams.dzU; obj.zUa = cParams.zUa;
            obj.lambda = cParams.lambda; obj.dlam = cParams.dlam; 
            obj.xL = cParams.xL; obj.sL = cParams.sL;
            obj.tau = cParams.tau;
            obj.xU = cParams.xU; obj.sU = cParams.sU;
            obj.bU = cParams.bU; obj.bL = cParams.bL;
            obj.gradient = cParams.gradient; obj.jacobian = cParams.jacobian; obj.residualLS = cParams.residualLS;
            obj.H = cParams.H;
            obj.filter = cParams.filter; obj.store = cParams.store; obj.status = cParams.status;
            obj.invDiagdU = cParams.invDiagdU; obj.invDiagdL = cParams.invDiagdL;
            obj.th = cParams.th; obj.ph = cParams.ph; obj.thetaMax = cParams.thetaMax;
            obj.sigmaL = cParams.sigmaL; obj.sigmaU = cParams.sigmaU;
            obj.iter = cParams.iteration;
            obj.firstLHS = cParams.firstLHS;
            obj.e = cParams.e;
            obj.lam = cParams.lam;
        end

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
    end

    methods (Static, Access = private)
        function meritFunction = meritFunctionComputer(cParams)
            merit = MeritComputer(cParams);
            merit.compute();
            meritFunction = merit.merit;
        end
    end
end