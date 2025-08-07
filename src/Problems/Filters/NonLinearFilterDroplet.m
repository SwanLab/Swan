classdef NonLinearFilterDroplet < handle
    
    properties (Access = private)
        mesh
        trial
        epsilon
        alpha
    end

    properties (Access = private)
        lambda
        direction
        sVar
        M
        K
        chiN
        proxdN
    end

    methods (Access = public)
        function obj = NonLinearFilterDroplet(cParams)
            obj.init(cParams);
            obj.createDirection(cParams);
            obj.updatePreviousGuess(0);
            obj.createMassMatrix();
            obj.createStiffnessMatrix();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            obj.createRHSChi(fun,quadOrder);
            iter = 1;
            tolerance = 1;
            while tolerance >= 1e-4 && iter<=25
                oldRho = obj.trial.fValues;
                obj.createRHSShapeDerivative(quadOrder);
                obj.solveProblem();
                obj.updatePreviousGuess(iter);
                tolerance = norm(obj.trial.fValues - oldRho)/norm(obj.trial.fValues);
                iter = iter + 1;
             end
             xF.setFValues(obj.trial.fValues);
        end

        function updateEpsilon(obj,eps)
            obj.alpha   = (obj.alpha/obj.epsilon)*eps;
            obj.lambda  = obj.lambda*(obj.epsilon/eps);
            obj.epsilon = eps;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, 1, 'P1'); % rho_eps
            obj.mesh    = cParams.mesh;
            obj.epsilon = obj.mesh.computeMeanCellSize();
            obj.alpha   = cParams.alpha*obj.epsilon;
            obj.lambda  = 10;
        end

        function createDirection(obj,s)
            th            = s.theta;
            k             = [cosd(th);sind(th)];
            obj.direction = ConstantFunction.create(k,obj.mesh);
        end

        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            LHS     = LHSIntegrator.create(s);
            obj.M   = LHS.compute();
        end

        function createStiffnessMatrix(obj)
            s.type  = 'StiffnessMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            LHS     = LHSIntegrator.create(s);
            obj.K   = LHS.compute();
        end

        function createRHSChi(obj,fun,quadType)
            switch class(fun)
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    s.mesh = fun.unfittedMesh;
                    s.type = 'Unfitted';
                otherwise
                    s.mesh = obj.mesh;
                    s.type = 'ShapeFunction';
            end
            s.quadType = quadType;
            int        = RHSIntegrator.create(s);
            test       = obj.trial;
            intN       = int.compute(fun,test);
            obj.chiN   = intN;
        end

        function createRHSShapeDerivative(obj,quadOrder)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            s.test     = obj.trial;
            int        = RHSIntegrator.create(s);
            rhs        = int.compute(obj.sVar);
            obj.proxdN = rhs;
        end

        function g = computeGradient(obj)
            s      = obj.sVar;
            a      = obj.alpha;
            maxFun = obj.CreateDomainMax(s);
            g      = s./(2*obj.lambda)-(a^2).*maxFun;
        end

        function m = CreateDomainMax(obj,sFun)
            s.operation = @(xV) max(zeros(size(xV(1,:,:))),sFun.evaluate(xV));
            m           = DomainFunction(s);
        end

        function solveProblem(obj)
            eps = obj.epsilon;
            l   = obj.lambda;
            LHS = obj.M + (eps^2/l)*obj.K;
            RHS = obj.chiN+(eps^2/l)*obj.proxdN;
            rhoi = LHS\RHS;
            obj.trial.setFValues(rhoi);
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
                s.mesh      = obj.mesh;
                A = DomainFunction(s);

                k        = obj.direction;
                gRhoK   = DP(gradRho,k);

                s.operation = @(xV) (squeezeParticular(gradRho.evaluate(xV),2)-A.evaluate(xV).*gRhoK.evaluate(xV).*k.evaluate(xV))./(1+2*l*muEv(xV));
                obj.sVar = DomainFunction(s);
            end
        end

        function mu = computeMu(obj)
            k       = obj.direction;
            l       = obj.lambda;
            a       = obj.alpha;
            gradRho = Grad(obj.trial);
            gRhoK   = DP(gradRho,k);

            % Mu = 0
            l2gRho  = Norm(gradRho,'L2');
            A       = a*sqrt(1+4*l+4*l^2*a^2)/(1+2*l*a^2);
            s.operation = @(xV) gRhoK.evaluate(xV)>=l2gRho/A;
            s.mesh     = obj.mesh;
            MuZeroCond = DomainFunction(s);

            % Mu = 1
            s.operation = @(xV) gRhoK.evaluate(xV)<=l2gRho/a;
            MuOneCond = DomainFunction(s);

            muGen = obj.createGeneralMu(l2gRho,gRhoK);

            s.operation = @(xV) 0.*MuZeroCond.evaluate(xV) + ...
                                1.*(MuOneCond.evaluate(xV) & not(MuZeroCond.evaluate(xV))) + ...
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
            s.mesh      = obj.mesh;
            mu          = DomainFunction(s);
        end
    end
end