% BPOPT Nonlinear Programming Solver
% Large-scale Interior Point solution methods
% Please note that this development version in MATLAB is primarily used
%   for algorithm development and small (<1000 variables) problems. An
%   optimized and large-scale version of this algorithm is available for
%   use at http://www.APMonitor.com and http://www.APOPT.com.
%
% Inputs
%   bp  = problem from bp_create
%
% Outputs
%   sol = solution returned as a structure
classdef bpopt < handle
    properties (Access = public)
        sol
        x 
        xL
        xU
        bL
        bU
        s 
        initResidual
        bp
        gradient
        jacobian
        lambda
        okJacobian
        okHessian
        okObjJacobian
        alpha_pr
        alpha_du
        thetaMax
        thetaMin
        filter
        store
        hessian
        e
        th
        ph
    end
    properties (Access = private)
        m
        n 
        ns
        sL 
        sU 
        tau
        zL
        zU
    end

    methods (Access = public)
        function obj = bpopt(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.checkArguments();
            obj.loadInitialVariables();
            obj.checkConsistentBounds();
            obj.computeVariableConstraintMultiplyers();
            obj.computeConstraintMultiplyers();
            obj.checkDerivatives();
            obj.firstNormInfeasibilities();
            obj.iterationInitializer();
            obj.interiorPointMehtodsSolver();
            obj.displaySolution();
        end

    end
    methods (Access = private)
        function init(obj,cParams)
            obj.bp = cParams.bp;
        end

        function checkArguments(obj)
            if (nargin==0),
                disp('Insufficient arguments, call bp_create to create a problem');
                obj.sol = [];
                return
            end
        end

        function loadInitialVariables(obj)
            % initial variables
            obj.loadVariables()
            [obj.x,obj.xL,obj.xU,obj.bL,obj.bU] = bp_x_init(obj.bp);
            % add slack for inequality constraints only
            k = 0;
            m = max(size(obj.bL));
            si = [];
            obj.s = [];
            obj.sL = [];
            obj.sU = [];
            for i = 1:m,
                if(obj.bU(i) > obj.bL(i)),
                    k = k + 1;
                    si(k) = i;
                    obj.s(k) = 0;
                    obj.sL(k) = obj.bL(i);
                    obj.sU(k) = obj.bU(i);
                end
            end

            % system size
            % n variables
            obj.n = max(size(obj.x));
            % ns slack variables
            obj.ns = max(size(obj.s));

            % initial residuals
            obj.computeInitialResidual();
            % m constraints (inequality or equality)
            obj.m = max(size(obj.initResidual));
            % ones vector
            obj.e = ones(obj.n + obj.ns,1);

            % initial equation slack variables
            k = 0;
            for i = 1:obj.m,
                % no slack variable when bU(i) == bL(i)
                if(obj.bU(i) > obj.bL(i)),
                    k = k + 1;
                    if (obj.bp.slack_init),
                    obj.s(k) = obj.initResidual(i);  % corrected, was -r(i)
                    else
                    obj.s(k) = 0.01;
                    end
                end
            end

            % initial parameters
            obj.alpha_pr = 1.0;
            obj.alpha_du = 1.0;
            
        end

        function loadVariables(obj)
            u.bp = obj.bp;
            initV = bp_x_init(u);
            initV.create();
            obj.x = initV.x0C;
            obj.xL = initV.xLC;
            obj.xU = initV.xUC;
            obj.bU = initV.bUC;
            obj.bL = initV.bLC;
        end

        function computeInitialResidual(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            obj.initResidual = residualComputer(u);
        end

        function checkConsistentBounds(obj)
            % check for consistent bounds
            if (min(obj.xU - obj.xL)<0),
                disp('bounds error (xU < xL)')
                obj.sol = [];
                return
            end
            if (min(bU-bL)<0),
                disp('bounds error (bU < bL)')
                obj.sol = [];
                return
            end
        end

        function moveVariablesToFeasible(obj);
            % move x into feasible region and off boundary initially
            for i = 1:obj.n,
                if (obj.x(i) <= obj.xL(i) || obj.x(i) >= obj.xU(i)),
                    fprintf(1,'Moving x0 to interior region\n')
                    break
                end
            end
            % Move x variables to be feasible
            for i = 1:obj.n,
                if (obj.x(i) <= obj.xL(i)),
                    obj.x(i) = min(obj.xU(i),obj.xL(i)+1e-2);
                end
                if (obj.x(i) >= obj.xU(i)),
                    obj.x(i) = max(obj.xL(i),obj.xU(i)-1e-2);
                end
            end
            % Move slack variables to be feasible
            for i = 1:obj.ns,
                if (obj.s(i) <= obj.sL(i)),
                    obj.s(i) = min(obj.sU(i),obj.sL(i)+1e-2);
                end
                if (obj.s(i) >= obj.sU(i)),
                    obj.s(i) = max(obj.sL(i),obj.sU(i)-1e-2);
                end
            end
        end

        function computeVariableConstraintMultiplyers(obj)
            obj.tau = min(obj.bp.tau_max,100*obj.bp.mu);
            % variable constraint multipliers
            % zL*(x-xL) = mu  =>  zL = mu / (x-xL)
            %zL(1:n) = 1;
            for i = 1:obj.n
                zLC(i) = obj.bp.mu / (obj.x(i) - obj.xL(i));
            end
            for i = 1:obj.ns
                zLC(obj.n + i) = obj.bp.mu / (obj.s(i) - obj.sL(i));
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

        function computeConstraintMultiplyers(obj)
            obj.computeObjectiveGradient();
            obj.computeJacobian();

            lam = pinv(full(obj.jacobian*obj.jacobian'))*obj.jacobian*(obj.zL'- obj.zU'- obj.gradient');

            obj.lambda = lam';
        end

        function computeObjectiveGradient(obj)
            u.bp = obj.bo;
            u.x = obj.x;
            u.s =obj.s;
            grad = bp_objgrad(u);
            grad.compute();
            obj.gradient = grad.objGradient;
        end

        function computeJacobian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.bL = obj.bL;
            u.bU = obj.bU;
            jac = bp_jac(u);
            jac.compute();
            obj.jacobian = jac.pd;
        end

        function computeHessian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.lam = obj.lambda;
            hes = bp_hes(u)
            hes.compute();
            obj.hessian = hes.hess;
        end

        function checkDerivatives(obj)
            if (obj.bp.prob == 0),
                okJacobian = true;
                okHessian = true;
            else
                % test only for hard coded problems
                obj.verifyJacobian();
                obj.verifyHessian();
                % test objective function gradients
                obj.verifyObjectiveJacobian();
            end
        end

        function verifyJacobian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            vJac = bp_verify_jac(u);
            vJac.compute();
            obj.okJacobian = vJac.okJacobian;
        end

        function verifyHessian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            u.lam = obj.lambda;
            vHess = bp_verify_hes(u);
            vHess.compute();
            obj.okHessian = vHess.okHessian;
        end

        function verifyObjectiveJacobian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            vObjJac = bp_verify_objgrad(u);
            vObjJac.compute();
            obj.okObjJacobian = vObjJac.okFirstDerivative;
        end

        function firstNormInfeasibilities(obj)
            % 1-norm of infeasibilities
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            th = bp_theta(u);
            th.compute();
            theta = th.theta;
            obj.thetaMax = 10^4  * max(1,theta);
            obj.thetaMin = 10^-4 * max(1,theta);
        end

        function iterationInitializer(obj)
            % initialize filter
            u.bp = obj.bp;
            u.x = obj.x;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = oj.bU;
            phi = phiComputer(u);
            obj.filter = [obj.th_max phi];

            % initialize iteration count
            u.iter = 0;

            % print iteration zero
            obj.iterationPrinter();

            obj.store = [u.iter u.x obj.alpha_pr obj.alpha_du obj.bp.mu];
        end

        function iterationPrinter(obj)
            u.lam = obj.lambda;
            u.zL = obj.zL;
            u.zU = obj.zU;
            u.alpha_pr =obj.alpha_pr;
            u.alpha_du = oobj.alpha_du;
            u.s = obj.s;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.bL = obj.bL;
            u.bU = obj.bU;
            pr = bp_iprint(u);
            pr.print();
        end
        function interiorPointMehtodsSolver(obj)
            obj.iterationProcess();
        end
        function displaySolution(obj)
            %% Display final solution
            disp('Status')
            disp(obj.status)
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

            if (obj.bp.prob>=1),
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

            % record solution for function return
            obj.sol.status = obj.status;
            obj.sol.x = obj.x;
            obj.sol.s = obj.s;
            obj.sol.lam = obj.lam;
            obj.sol.zL = obj.zL;
            obj.sol.zU = obj.zU;
        end
    end

    methods (Static, Access = public)
        function residual = residualComputer(cParams);
            residualC = bp_res(cParams);
            residualC.compute();
            residual = res.c;
        end

        function phi = phiComputer(cParams)
            phiC = bp_phi(cParams);
            phiC.compute();
            phi = phiC.phi;
            obj.ph = phi;
        end

        function theta = thetaComputer(cParams)
            th = bp_theta(cParams);
            th.compute();
            theta = th.theta;
            obj.th = theta;
        end
    end
end