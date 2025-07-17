classdef NonLinearFilterSegment < handle
    
    properties (Access = private)
        mesh
        trial
        theta
        epsilon
        alpha
        beta
    end

    properties (Access = private)
        lineSearch
        direction
        sVar
        M
        K
        intChi
        rhsDer
    end

    methods (Access = public)
        function obj = NonLinearFilterSegment(cParams)
            obj.init(cParams);
            obj.createDirection();
            obj.createMassMatrix();
            obj.createDirectionalStiffnessMatrix();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = copy(obj.trial);     
            obj.computeInitialGuess(fun,quadOrder);
            obj.intChi = obj.createRHSShapeFunction(fun,quadOrder);
            iter   = 1;
            error0 = inf;
            error  = inf;
            obj.updateDotProduct(obj.trial); 
            while error >= 1e-6  && iter<=1000
                obj.createRHSDirectionalDerivative(quadOrder);
                obj.solveProblem();
                obj.updateDotProduct(obj.trial);
                dJ = obj.computeCostGradient(quadOrder);
                error = Norm(dJ,'L2');
                iter = iter + 1;
                obj.updateLineSearch(error0,error);
                error0 = error;
            end
            % Projection ub, lb
            xF.setFValues(obj.trial.fValues);
        end

        function updateEpsilon(obj,eps)
            obj.alpha      = (obj.alpha/obj.epsilon)*eps;
            obj.beta       = (obj.beta/obj.epsilon)*eps;
            obj.lineSearch = obj.lineSearch*(obj.epsilon/eps);
            obj.epsilon    = eps;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial      = LagrangianFunction.create(cParams.mesh, 1, 'P1');
            obj.mesh       = cParams.mesh;
            obj.theta      = cParams.theta;
            obj.epsilon    = obj.mesh.computeMeanCellSize();
            obj.alpha      = cParams.alpha*obj.epsilon;
            obj.beta       = cParams.beta*obj.epsilon;
            obj.lineSearch = 10;
        end

        function createDirection(obj)
            th = obj.theta;            
            k  = [cosd(th);sind(th)];
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

        function createDirectionalStiffnessMatrix(obj)
            k       = obj.direction.constant;
            s.type  = 'AnisotropicStiffnessMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            s.A     = ConstantFunction.create(k*k',obj.mesh);
            LHS     = LHSIntegrator.create(s);
            obj.K   = LHS.compute();
        end

        function computeInitialGuess(obj,fun,quadOrder)
            rhoE1 = obj.trial;
            rhoE2 = fun.project(obj.trial.order);
            J1    = obj.computeCostFunction(rhoE1,fun,quadOrder);
            J2    = obj.computeCostFunction(rhoE2,fun,quadOrder);
            if J2<J1
                obj.trial.setFValues(rhoE2.fValues);
            end
        end

        function J = computeCostFunction(obj,rhoE,fun,quadOrder)
            obj.updateDotProduct(rhoE);
            h    = obj.computeMeasure();
            int1 = 0.5.*rhoE.^2 + 0.5.*(h.^2);
            int2 = fun.*(-rhoE);
            int3 = fun.*fun.*0.5;
            j1   = Integrator.compute(int1,obj.mesh,quadOrder);
            j2   = Integrator.compute(int2,obj.mesh,quadOrder);
            j3   = Integrator.compute(int3,obj.mesh,quadOrder);
            J    = j1 + j2 + j3;
        end

        function updateDotProduct(obj,rhoE)
            gradRho  = Grad(rhoE);
            k        = obj.direction;
            obj.sVar = DP(gradRho,k);
        end

        function h = computeMeasure(obj)
            s = obj.sVar;
            a = obj.alpha;
            b = obj.beta;
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            h = a.*maxFun - b.*minFun;
        end

        function createRHSDirectionalDerivative(obj,quadOrder)
            g          = obj.computeMeasureGradient();
            f          = obj.direction;
            obj.rhsDer = obj.createRHSShapeDerivative(f.*g,quadOrder);
        end

        function g = computeMeasureGradient(obj)
            s      = obj.sVar;
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            g      = s./(2*obj.lineSearch)-(obj.alpha^2).*maxFun-(obj.beta^2).*minFun;
        end

        function solveProblem(obj)
            tau  = obj.lineSearch;
            LHS  = obj.M + (1/(2*tau))*obj.K;
            RHS  = obj.rhsDer + obj.intChi;
            rhoi = LHS\RHS;
            obj.trial.setFValues(rhoi);
        end

        function dJ = computeCostGradient(obj,quadOrder)
            [rhs1,rhs2,rhs3] = obj.computeCostGateaux(quadOrder);
            Ms   = obj.M;
            Ms   = diag(sum(Ms,1));
            fVal = Ms\(rhs1+rhs2+rhs3);
            dJ   = copy(obj.trial);
            dJ.setFValues(fVal);
        end

        function [rhs1,rhs2,rhs3] = computeCostGateaux(obj,quadOrder)
            a      = obj.alpha;
            b      = obj.beta;
            s      = obj.sVar;
            k      = obj.direction;
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            rhs1   = obj.createRHSShapeFunction(obj.trial,quadOrder);
            rhs2   = -obj.intChi;
            int3   = (a^2.*maxFun + b^2.*minFun).*k;
            rhs3   = obj.createRHSShapeDerivative(int3,quadOrder);
        end

        function updateLineSearch(obj,error0,error)
            if error<=error0
                obj.lineSearch = obj.lineSearch*1.1;
            else
                obj.lineSearch = obj.lineSearch/2;
            end
        end

        function intN = createRHSShapeFunction(obj,fun,quadType)
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
        end

        function intdN = createRHSShapeDerivative(obj,fun,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSIntegrator.create(s);
            test       = obj.trial;
            intdN      = int.compute(fun, test);
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
    end
end