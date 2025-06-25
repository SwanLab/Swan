classdef NonLinearFilterSegment < handle
    
    properties (Access = private)
        mesh
        trial
        sVar
        theta
        alpha
        beta
    end

    properties (Access = private)
        lineSearch
        direction
        directionFunction
        M
        K
        intChi
        rhsDer
        rhoDif
        Den
        a2
        b2
        rhoOld
        epsilonOld
    end

    methods (Access = public)

        function obj = NonLinearFilterSegment(cParams)
            obj.init(cParams);
            obj.createDirection();
            obj.createDirectionFunction();
            obj.updateDotProductPreviousGuess();
            obj.createMassMatrix();
            obj.createDirectionalStiffnessMatrix();
            obj.Den = obj.createConstantFunction(2*obj.lineSearch);
            eps = obj.mesh.computeMeanCellSize();
            obj.epsilonOld = eps;
            obj.a2  = obj.createConstantFunction((obj.alpha*eps)^2);
            obj.b2  = obj.createConstantFunction((obj.beta*eps)^2);
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);     
            obj.trial = fun.project(obj.rhoOld.order);
            obj.createRHSChi(fun,quadOrder);
            iter = 1;
            tolerance = 1;
            obj.updateDotProductPreviousGuess();      
            while tolerance >= 1e-4  && iter<=1000
                obj.rhoOld.setFValues(obj.trial.fValues);
                obj.createRHSDirectionalDerivative(quadOrder);
                obj.solveProblem();
                obj.updateDotProductPreviousGuess();
                obj.rhoDif.setFValues(obj.rhoOld.fValues - obj.trial.fValues);
                tolerance = Norm(obj.rhoDif,'L2')/Norm(obj.trial,'L2');
                iter = iter + 1;
            end
            xF.setFValues(obj.trial.fValues);
        end

        function updateEpsilon(obj,eps)
            obj.a2  = obj.createConstantFunction((obj.alpha*eps)^2);
            obj.b2  = obj.createConstantFunction((obj.beta*eps)^2);
            obj.lineSearch = obj.lineSearch*(obj.epsilonOld/eps);
            obj.Den = obj.createConstantFunction(2*obj.lineSearch);
            obj.epsilonOld = eps;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial = LagrangianFunction.create(cParams.mesh, 1, 'P1');
            obj.rhoOld = LagrangianFunction.create(cParams.mesh, 1, 'P1');
            obj.rhoDif = LagrangianFunction.create(cParams.mesh, 1, 'P1');
            obj.mesh  = cParams.mesh;
            obj.theta = cParams.theta;
            obj.alpha = cParams.alpha;
            obj.beta  = cParams.beta;
            obj.lineSearch = 2e-2;
        end

        function createDirection(obj)
            th = obj.theta;            
            k = [cosd(th);sind(th)];
            obj.direction = k;
        end

        function createDirectionFunction(obj)
            k = obj.direction;
            obj.directionFunction  = obj.createConstantFunction(k);
        end

        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            LHS     = LHSIntegrator.create(s);
            obj.M   = LHS.compute();
        end

        function createDirectionalStiffnessMatrix(obj)
            k       = obj.direction;
            s.type  = 'AnisotropicStiffnessMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            s.aniAlphaDeg = 0;
            s.CAnisotropic = k*k';
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
            rhs        = int.compute(fun,test);
            obj.intChi = rhs;
        end

        function createRHSDirectionalDerivative(obj,quadOrder)
            s.mesh = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSIntegrator.create(s);
            test       = obj.trial;
            g          = obj.computeGradient();
            f          = obj.directionFunction;
            rhs        = int.compute(f.*g, test);
            obj.rhsDer = rhs;
        end


        function g = computeGradient(obj)
            s   = obj.sVar;
%             Den = obj.createConstantFunction(2*obj.lineSearch);
%             a2  = obj.createConstantFunction(a^2);
%             b2  = obj.createConstantFunction(b^2);
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            g   = s./obj.Den-obj.a2.*maxFun-obj.b2.*minFun;
        end

        function m = createDomainMax(obj,sFun)
            s.operation = @(xV) max(zeros(size(xV(1,:,:))),sFun.evaluate(xV));
            s.mesh      = obj.mesh;
            m           = DomainFunction(s);
        end

        function m = createDomainMin(obj,sFun)
            s.operation = @(xV) min(zeros(size(xV(1,:,:))),sFun.evaluate(xV));
            s.mesh      = obj.mesh;
            m           = DomainFunction(s);
        end

        function solveProblem(obj)
            tau = obj.lineSearch;
            LHS = obj.M + (1/(2*tau))*obj.K;
            RHS = obj.rhsDer + obj.intChi;
            rhoi = LHS\RHS;
            obj.trial.setFValues(rhoi);
        end

        function aF = createConstantFunction(obj,c)
            s.ndimf = length(c);
            s.operation = @(xV) obj.createOperation(c,xV);
            s.mesh = obj.mesh;
            aF = DomainFunction(s);
        end

        function b = createOperation(obj,c,xV)
            b = c.*ones([length(c),size(xV,2),obj.mesh.nelem]); 
        end

        function updateDotProductPreviousGuess(obj)
            gradRho  = Grad(obj.trial);
            k        = obj.directionFunction;
            obj.sVar = DP(gradRho,k);
%             obj.sVar = Integrator.compute(gradRho.*k,obj.mesh,2);
        end
    end
end