classdef LineSearchComputer < InteriorPointMethodsSolver

    properties (Access = public)
        e_mu
    end
    properties (Access = private)
        predictReduction
    end
    methods (Access = public)
        function search(obj)
            obj.createLineSearchMethod();
            obj.updateVariables();
            obj.updateZvalues();
            obj.pushVariables();
            obj.reduceMeritFunction();
            obj.performLineSearch();
            obj.checkAcceptance();
            obj.updateWithAcceptance();
            obj.checkConvergence();
            obj.checkTerminationConditions();
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
                    obj.alpha_du = alpha_du_max;
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
                obj.lambda = obj.lambda + obj.alpha_pr * obj.dlam';
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
            obj.predictReduction = -obj.alpha_pr*obj.gradient*[obj.dx;obj.ds] - 0.5*obj.alpha_pr^2*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds] + obj.bp.nu*(norm(obj.residual',1) - norm(obj.residual' + obj.alpha_pr*obj.jacobian*[obj.dx;obj.ds],1));
            u.bp = obj.bp;
            u.x = obj.x;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.s = obj.s;
            u.bU = obj.bU;
            u.bL = obj.bL;
            meritX = meritFunctionComputer(u);
            
            u.x = obj.xa;
            u.s = obj.sa;
            meritA = meritFunctionComputer(u);
            obj.reduced = meritX - meritA;
        end

        function performLineSearch(obj)
            % line search criteria
            %  1 = reduction in merit function
            %  2 = simple clipping
            %  3 = filter method
            switch obj.bp.line_search
            case(1)
                % merit function
                pred_phi_1 = obj.gradient*[obj.dx;obj.ds];
                % set 2nd derivative contribution to zero if < 0
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
                %if (new_nu>bp.nu),
                %    bp.nu = min(1000,new_nu + 1);
                %end
                obj.bp.nu = max(1,min(1000,nuUpdated));
                % update predicted and actual reductions with new nu value
                obj.predictUpdated = -obj.alpha_pr*obj.gradient*[obj.dx;obj.ds] - 0.5*obj.alpha_pr^2*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds] + obj.bp.nu*(norm(obj.residual',1) - norm(obj.residual' + obj.alpha_pr*obj.jacobian*[obj.dx;obj.ds],1));
                u.bp = obj.bp;
                u.x = obj.x;
                u.xL = obj.xL;
                u.xU = obj.xU;
                u.s = obj.s;
                u.bU = obj.bU;
                u.bL = obj.bL;
                meritX = meritFunctionComputer(u);
                
                u.x = obj.xa;
                u.s = obj.sa;
                meritA = meritFunctionComputer(u);
                obj.reducedUpdated = meritX - meritA;
                eta = 0.2;
                % compare actual reduction to predicted reduction
                % as long as the actual reduction is a fraction of the
                % predicted reduction then accept the trial point
                if (obj.reducedUpdated >= eta*obj.predictUpdated)
                    obj.accept = true;
                else
                    obj.accept = false;
                end
            case(2)
                obj.accept = true;
                % test case for alpha_pr, alpha_du = 1 with clipping only
                obj.alpha_pr = 1;
                obj.alpha_du = 1;
            case(3)
                % check if acceptable point with filter method
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
        function checkAcceptance(obj)
            % testing for acceptance
            %ac = true;
            % second order correction
            if (false)
                if (not(obj.accept))
                    % new residuals
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
                    % compute search direction, solving Ax = b
                    d3 = -pinv(full(obj.firstLHS))*b1;
                    dx3 = d3(1:obj.n,1);
                    ds3 = d3(obj.n+1:obj.n+obj.ns,1);
                    dlam3 = d3(obj.n + obj.ns + 1:obj.n + obj.ns + obj.m,1);
                    % compute search direction for z (explicit solution)
                    dzL3 = obj.bp.mu*obj.invDiagdL*obj.e - obj.zL' - obj.sigmaL*[dx3; ds3];
                    dzU3 = obj.bp.mu*obj.invDiagdU*obj.e - obj.zU' + obj.sigmaU*[dx3; ds3];
                        
                    % compute new acceptance point
                    obj.xa = obj.x + obj.alpha_pr * dx3';
                    obj.sa = obj.s + obj.alpha_pr * ds3';
                    obj.lambda = obj.lambda + obj.alpha_du * dlam3';
                    obj.zLa = obj.zL + obj.alpha_du * dzL3';
                    obj.zUa = obj.zU + obj.alpha_du * dzU3';
                        
                    % check for acceptance of trial point
                    u.bp = obj.bp;
                    u.x = obj.x;
                    u.xa = obj.xa;
                    u.xL = obj.xL;
                    u.xU = obj.xU;
                    u.filter = obj.filter;
                    
                    ac = bp_accept(u);
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
                % accept point
                obj.x = obj.xa;
                obj.s = obj.sa;
                obj.lam = obj.lama;
                obj.zL = obj.zLa;
                obj.zU = obj.zUa;
                
                % update filter
                obj.filter = [obj.filter; obj.th obj.ph];
            end
        end
        function checkConvergence(obj)
                % check for convergence
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
                residual = residualComputer(u);
                part(2) = max(abs(residual));
                % use mu = 0 to test for convergence
                part(3) = max(abs(diag([obj.x-obj.xL obj.s-obj.sL])*diag(obj.zL)*obj.e))/s_c;
                part(4) = max(abs(diag([obj.xU-obj.x obj.sU-obj.s])*diag(boj.zU)*obj.e))/s_c;
                obj.e_mu = max(part);
                        
                if (obj.bp.idebug>=1)
                    fprintf(1,'pred: %12.4e, ared: %12.4e, e_mu: %12.4e, k*mu: %12.4e\n',obj.predictUpdated,obj.reducedUpdated,obj.e_mu,0.2*obj.bp.mu);
                end
        end
        function checkTerminationConditions(obj)
            % check for termination conditions
            if (obj.e_mu <= obj.bp.e_tol)
                fprintf(1,'\nSuccessful solution\n');
                obj.status = 'success';
                break;
            end
        end
        function checkNewBarrierProblem(obj)
            % check for new barrier problem
            k_mu  = 0.2; % (0,1)
                
            % only update mu if requested (not debugging)
            if (obj.bp.mu_update)
                if (obj.e_mu < k_mu * obj.bp.mu)
                    th_mu = 1.5; % (1,2)
                    % update mu
                    obj.bp.mu = max(bp.e_tol/10,min(k_mu*obj.bp.mu,obj.bp.mu^th_mu));
                    % update tau
                    tau = min(obj.bp.tau_max,100*obj.bp.mu);
                    % re-initialize filter
                    u.bp = obj.bp;
                    u.x = obj.x;
                    u.xL = obj.xL;
                    u.s = obj.s;
                    u.bL = obj.bL;
                    u.bU = obj.bU;
                    phi = phiComputer(u);
                    obj.filter = [obj.th_max phi];
                end
            end
        end
        function printAndStoreVariables(obj)
            % print iteration
            obj.printIteration();   
            obj.store = [obj.store; obj.iter obj.x obj.alpha_pr obj.alpha_du obj.bp.mu];
        end
        function printIteration(obj)
            u.bp = obj.bp;
            u.iter = obj.iter;
            u.x = obj.x;
            u.lambda = obj.lambda;
            u.zL = obj.zL;
            u.zU = obj.zU;
            u.alpha_pr = obj.alpha_pr;
            u.alpha_du = obj.alpha_du;
            u.s = obj.s;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.xL = obj.xL;
            u.bL = obj.bL;
            u.bU = obj.bU;
            bp_iprint(u);
        end 
        function checkMaxIter(obj)
            % reached the end of the iteration loop without convergence
            if (obj.iter == obj.bp.maxiter)
                obj.status = 'failed: max iterations';
                break;
            end
            
           
        end
        function checkNotAcceptance(obj)
            % don't do a line search in this MATLAB version
            % just cycle through on another iteration with a lower alpha if
            %   the point was not accepted
            if (obj.accept)
                % reset alpha
                obj.alpha_pr = 1.0;
                obj.alpha_du = 1.0;
            else
                % reject point and move alpha_x
                obj.alpha_pr = obj.alpha_pr / 2;
                obj.alpha_du = obj.alpha_du / 2;
                if (obj.alpha_pr < 1e-4)
                    obj.alpha_pr = 1.0;
                    obj.alpha_du = 1.0;
                end
            end
        end
    end

    methods (Static, Access = public)
        function meritFunction = meritFunctionComputer(cParams)
            merit = bp_merit(cParams);
            merit.compute();
            meritFunction = merit.merit;
        end
    end
end