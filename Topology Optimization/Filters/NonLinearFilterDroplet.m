classdef NonLinearFilterDroplet < handle
    
    properties (Access = private)
        mesh
        epsilon
        trial
        sVar
        theta
        alpha
        beta
    end

    properties (Access = private)
        lambda
        direction
        M
        K
        intChi
        rhsDer
        Den
        a2
        b2
    end

    methods (Access = public)
        function obj = NonLinearFilterDroplet(cParams)
            obj.init(cParams);
            obj.createDirection();
            obj.updatePreviousGuess(0);
            obj.createMassMatrix();
            obj.createStiffnessMatrix();
%             obj.Den = obj.createConstantFunction(2*obj.lambda);
%             obj.a2  = obj.createConstantFunction(obj.alpha^2);
%             obj.b2  = obj.createConstantFunction(obj.beta^2);
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            obj.createRHSChi(fun,quadOrder);
            iter = 1;
            tolerance = 1;
%             filename = 'beta5alpha5theta30.gif';
            while tolerance >= 1e-5 
                oldRho = obj.trial.fValues;
                obj.createRHSShapeDerivative(quadOrder);
                obj.solveProblem();
                obj.updatePreviousGuess(iter);
                tolerance = norm(obj.trial.fValues - oldRho)/norm(obj.trial.fValues);
                iter = iter + 1;
                 %disp(iter);  
                 disp(tolerance);
             end
           
            obj.trial.plot
            xF.fValues = obj.trial.fValues;

        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial = LagrangianFunction.create(cParams.mesh, 1, 'P1'); % rho_eps
            obj.mesh  = cParams.mesh;
            obj.epsilon = cParams.epsilon;
            obj.theta = cParams.theta;
            obj.alpha = cParams.alpha;
            obj.beta  = 0;
            obj.lambda = 1;
        end

        function createDirection(obj)
            th            = obj.theta;
            obj.direction = [cosd(th);sind(th)];
        end

        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
        end

        function createStiffnessMatrix(obj)
            s.type  = 'StiffnessMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            LHS     = LHSintegrator.create(s);
            obj.K   = LHS.compute();
        end

        function createRHSChi(obj,fun,quadOrder)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(fun,test);
            obj.intChi = rhs;
        end

        function createRHSShapeDerivative(obj,quadOrder)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(obj.sVar, test);
            obj.rhsDer = rhs;
        end

        function aF = createAnalyticalDirection(obj)
            k = obj.direction;
            s.fHandle = @(x) [k(1)*ones(size(x(1,:,:)));k(2)*ones(size(x(1,:,:)))];
            s.ndimf = 2;
            s.mesh = obj.mesh;
            aF = AnalyticalFunction(s);
        end

        function g = computeGradient(obj)
            s   = obj.sVar;
%             Den = obj.createConstantFunction(2*obj.lineSearch);
%             a2  = obj.createConstantFunction(a^2);
%             b2  = obj.createConstantFunction(b^2);
            maxFun = obj.CreateDomainMax(s);
            minFun = obj.CreateDomainMin(s);
            g   = s./obj.Den-obj.a2.*maxFun-obj.b2.*minFun;
        end

        function m = CreateDomainMax(obj,sFun)
            s.operation = @(xV) max(zeros(size(xV(1,:,:))),sFun.evaluate(xV));
            m           = DomainFunction(s);
        end

        function m = CreateDomainMin(obj,sFun)
            s.operation = @(xV) min(zeros(size(xV(1,:,:))),sFun.evaluate(xV));
            m           = DomainFunction(s);
        end

        function solveProblem(obj)
            eps = obj.epsilon;
            l   = obj.lambda;
            LHS = obj.M + (eps^2/l)*obj.K;
            RHS = obj.intChi+(eps^2/l)*obj.rhsDer;
            rhoi = LHS\RHS;
            obj.trial.fValues = rhoi;
        end

        function aF = createConstantFunction(obj,c)
            s.fHandle = @(x) c*ones(size(x(1,:,:)));
            s.ndimf = 1;
            s.mesh = obj.mesh;
            aF = AnalyticalFunction(s);
        end

        function updatePreviousGuess(obj,iter)
            if iter == 0
                obj.sVar = LagrangianFunction.create(obj.mesh, 2, 'P0');
            else
                l        = obj.lambda;
                a        = obj.alpha;
                mu       = obj.computeMu();
                gradRho  = Grad(obj.trial);

                muEv = @(xV) mu.evaluate(xV);
                s.operation = @(xV) 2*l*(1-muEv(xV))*a^2./(1+2*l*muEv(xV)+2*l*(1-muEv(xV))*a^2);
                A = DomainFunction(s);

                k        = obj.createAnalyticalDirection();
                gRhoK   = DotProduct(gradRho,k);

                s.operation = @(xV) (gradRho.evaluate(xV)-A.evaluate(xV).*gRhoK.evaluate(xV).*k.evaluate(xV))./(1+2*l*muEv(xV));
                obj.sVar = DomainFunction(s);
            end
        end

        function mu = computeMu(obj)
            k       = obj.createAnalyticalDirection();
            l       = obj.lambda;
            a       = obj.alpha;
            gradRho = Grad(obj.trial);
            gRhoK   = DotProduct(gradRho,k);

            % Mu = 0
            l2gRho  = L2Norm.compute(obj.mesh,gradRho);
            A       = a*sqrt(1+4*l+4*l^2*a^2)/(1+2*l*a^2);
            s.operation = @(xV) gRhoK.evaluate(xV)>=l2gRho/A;
            MuZeroCond = DomainFunction(s);

            % Mu = 1
            s.operation = @(xV) gRhoK.evaluate(xV)<=l2gRho/a;
            MuOneCond = DomainFunction(s);

            muGen = obj.createGeneralMu(l2gRho,gRhoK);

            s.operation = @(xV) 0.*MuZeroCond.evaluate(xV) + ...
                                1.*MuOneCond.evaluate(xV) + ...
                                muGen.evaluate(xV).*(not(MuZeroCond.evaluate(xV)) & not(MuOneCond.evaluate(xV)));
            mu = DomainFunction(s);
        end

        function mu = createGeneralMu(obj,l2gRho,gRhoK)
            l     = obj.lambda;
            a     = obj.alpha;
            gRKEv = @(xV) gRhoK.evaluate(xV);

            A = @(xV) (1+2*l*a^2)*sqrt(l2gRho^2-(gRKEv(xV)).^2);
            B = @(xV) sqrt(a^2-1)*gRKEv(xV);
            C = 2*l*sqrt(a^2-1);
            D = @(xV) gRKEv(xV)+sqrt((a^2-1)*(l2gRho^2-(gRKEv(xV)).^2));
            D = @(xV) D(xV).*(D(xV)~=0)+1.*(D(xV)==0);

            s.operation = @(xV) (A(xV)-B(xV))./(C*D(xV));
            mu          = DomainFunction(s);
        end
    end
end