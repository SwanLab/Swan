classdef OptimizerComputer < InteriorPointMethodsSolver
    properties (Access = public)
    end
    properties (Access = protected)
    end

    methods (Access = public)
        
        function compute(obj)
            obj.checkArguments();
            obj.loadInitialVariables();
            obj.checkConsistentBounds();
            obj.computeVariableConstraintMultiplyers();
            obj.computeConstraintMultiplyers();
            obj.checkDerivatives();
            obj.firstNormInfeasibilities();
            obj.iterationInitializer();
            obj.iterationProcess();
            %obj.displaySolution();
        end

    end
    methods (Access = protected)
        function init(obj,cParams)
            obj.bp = cParams;
        end

        function checkArguments(obj)
            if (nargin==0)
                disp('Insufficient arguments, call bp_create to create a problem');
                obj.sol = [];
                return
            end
        end

        function loadInitialVariables(obj)
            % initial variables
            obj.loadVariables()
            % add slack for inequality constraints only
            k = 0;
            obj.m = max(size(obj.bL));
            si = [];
            obj.s = [];
            obj.sL = [];
            obj.sU = [];
            for i = 1:obj.m
                if(obj.bU(i) > obj.bL(i))
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
            for i = 1:obj.m
                % no slack variable when bU(i) == bL(i)
                if(obj.bU(i) > obj.bL(i))
                    k = k + 1;
                    if (obj.bp.slack_init)
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
            u = obj.bp;
            initV = ParametersInitializer(u);
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
            obj.initResidual = obj.residualComputer(u);
        end

        function checkConsistentBounds(obj)
            % check for consistent bounds
            if (min(obj.xU - obj.xL)<0)
                disp('bounds error (xU < xL)')
                obj.sol = [];
                return
            end
            if (min(obj.bU - obj.bL)<0)
                disp('bounds error (bU < bL)')
                obj.sol = [];
                return
            end
        end

        function moveVariablesToFeasible(obj)
            % move x into feasible region and off boundary initially
            for i = 1:obj.n
                if (obj.x(i) <= obj.xL(i) || obj.x(i) >= obj.xU(i))
                    fprintf(1,'Moving x0 to interior region\n')
                    break
                end
            end
            % Move x variables to be feasible
            for i = 1:obj.n
                if (obj.x(i) <= obj.xL(i))
                    obj.x(i) = min(obj.xU(i),obj.xL(i)+1e-2);
                end
                if (obj.x(i) >= obj.xU(i))
                    obj.x(i) = max(obj.xL(i),obj.xU(i)-1e-2);
                end
            end
            % Move slack variables to be feasible
            for i = 1:obj.ns
                if (obj.s(i) <= obj.sL(i))
                    obj.s(i) = min(obj.sU(i),obj.sL(i)+1e-2);
                end
                if (obj.s(i) >= obj.sU(i))
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

        function checkDerivatives(obj)
            if (obj.bp.prob == 0)
                obj.okJacobian = true;
                obj.okHessian = true;
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
            vJac = JacobianVerifier(u);
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
            vHess = HessianVerifier(u);
            vHess.compute();
            obj.okHessian = vHess.okHessian;
        end

        function verifyObjectiveJacobian(obj)
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            vObjJac = GradientVerifier(u);
            vObjJac.check();
            obj.okObjJacobian = vObjJac.okFirstDerivative;
        end

        function firstNormInfeasibilities(obj)
            % 1-norm of infeasibilities
            u.bp = obj.bp;
            u.x = obj.x;
            u.s = obj.s;
            u.bL = obj.bL;
            u.bU = obj.bU;
            the = ThetaComputer(u);
            the.compute();
            theta = the.theta;
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
            u.bU = obj.bU;
            phi = obj.phiComputer(u);
            obj.filter = [obj.thetaMax phi];
            % initialize iteration count
            iter = 0;
            % print iteration zero
            obj.iterationPrinter();
            obj.store = [iter u.x obj.alpha_pr obj.alpha_du obj.bp.mu];
        end

        function iterationPrinter(obj)
            u.lam = obj.lambda;
            u.zL = obj.zL;
            u.zU = obj.zU;
            u.alpha_pr =obj.alpha_pr;
            u.alpha_du = obj.alpha_du;
            u.s = obj.s;
            u.x = obj.x;
            u.xL = obj.xL;
            u.xU = obj.xU;
            u.bL = obj.bL;
            u.bU = obj.bU;
            u.bp = obj.bp;
            u.iter = 0;
            pr = PrinterComputer(u);
            pr.print();
        end
    end
end