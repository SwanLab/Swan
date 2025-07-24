classdef NonLinearFilterSegment < handle

    properties (Access = public)
            errorVec
            iterVec
    end

    properties (Access = private)
        mesh
        trial
        theta
        epsilon
        alpha
        beta
        ub
        lb
    end

    properties (Access = private)
        lineSearch
        filter
        direction
        sVar % Eliminable?
        M
        K % Mergeable?
        intChi
        rhsDer % Mergeable?
        isBoundFree
    end

    methods (Access = public)
        function obj = NonLinearFilterSegment(cParams)
            obj.init(cParams);
            obj.createFilterInitialGuess();
            obj.createDirection();
            obj.createMassMatrix();
            obj.createDirectionalStiffnessMatrix();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = copy(obj.trial);     
            obj.computeInitialGuess(fun,quadOrder);
            obj.updateDotProduct(obj.trial); 
            obj.intChi = obj.createRHSShapeFunction(fun,quadOrder);
            iter   = 1;
            dJ0 = obj.computeCostGradient(quadOrder);
            error0 = Norm(dJ0,'L2');
            error  = inf;
            while error >= 1e-6  && iter<=1000
                x0 = copy(obj.trial);
                isAcceptable = false;
                while not(isAcceptable)
                    obj.createRHSDirectionalDerivative(quadOrder);
                    obj.solveProblem();
                    obj.updateDotProduct(obj.trial);
                    dJ = obj.computeCostGradient(quadOrder);
                    error = Norm(dJ,'L2');
                    if error <= error0
                        obj.lineSearch = obj.lineSearch*1.1;
                        isAcceptable = true;
                    else
                        obj.lineSearch = obj.lineSearch/2;
                        obj.trial.setFValues(x0.fValues);
                        obj.updateDotProduct(obj.trial);
                    end
                end
                obj.iterVec = [obj.iterVec; iter];
                obj.errorVec = [obj.errorVec; error];
                iter = iter + 1;
                error0 = error;
            end
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
            obj.ub         = cParams.ub;
            obj.lb         = cParams.lb;
            obj.lineSearch = 10;
        end

        function createFilterInitialGuess(obj)
            s.filterType = 'LUMP';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,obj.trial.order);
            f = Filter.create(s);
            obj.filter = f;
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
            rhoE2 = obj.filter.compute(fun,quadOrder);
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
            e = obj.mesh.computeMeanCellSize();
            s = obj.sVar;
            a = max(obj.alpha,0);
            b = max(obj.beta,0);
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
            e      = obj.mesh.computeMeanCellSize();
            s      = obj.sVar;
            a      = max(obj.alpha,0);
            b      = max(obj.beta,0);
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            g      = s./(2*obj.lineSearch)-(a^2).*maxFun-(b^2).*minFun;
        end

        function solveProblem(obj)
            tau  = obj.lineSearch;
            LHS  = obj.M + (1/(2*tau))*obj.K;
            RHS  = obj.rhsDer + obj.intChi;
            x    = LHS\RHS;
            %rhoi = max(obj.lb,min(obj.ub,x));
            rhoi = x;
            obj.trial.setFValues(rhoi);
            obj.isBoundFree = find((x-rhoi)==0);
        end

        function dJ = computeCostGradient(obj,quadOrder)
            dJUnc = obj.computeUnconstrainedCostGradient(quadOrder);
            dj    = copy(dJUnc);
            djV   = dj.fValues;
            djV(obj.isBoundFree) = 0;
            dj.setFValues(djV);
            mu = obj.createDomainMax(-dj);
            l  = obj.createDomainMax(dj);
            %dJ = dJUnc + mu - l;
            dJ = dJUnc;
        end

        function dJ = computeUnconstrainedCostGradient(obj,quadOrder)
            [rhs1,rhs2,rhs3] = obj.computeCostGateaux(quadOrder);
            Ms   = obj.M;
            Ms   = diag(sum(Ms,1));
            fVal = Ms\(rhs1+rhs2+rhs3);
            dJ   = copy(obj.trial);
            dJ.setFValues(fVal);
        end

        function [rhs1,rhs2,rhs3] = computeCostGateaux(obj,quadOrder)
            e      = obj.mesh.computeMeanCellSize();
            a      = max(obj.alpha,0);
            b      = max(obj.beta,0);
            s      = obj.sVar;
            k      = obj.direction;
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            rhs1   = obj.createRHSShapeFunction(obj.trial,quadOrder);
            rhs2   = -obj.intChi;
            int3   = (a^2.*maxFun + b^2.*minFun).*k;
            rhs3   = obj.createRHSShapeDerivative(int3,quadOrder);
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