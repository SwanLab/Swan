classdef NonLinearFilterSegment < handle

    properties (Access = private)
        mesh
        trial
        epsilon
        alphaEps
        betaEps
        alpha
        beta
        %ub
        %lb
    end

    properties (Access = private)
        lineSearch
        filter
        direction
        sVar
        M
        K
        chiN
        kdhdN
        %isBoundFree
        dNIntegrator
    end

    methods (Access = public)
        function obj = NonLinearFilterSegment(cParams)
            obj.init(cParams);
            obj.createFilterInitialGuess();
            obj.createDirection(cParams);
            obj.createMassMatrix();
            obj.createDirectionalStiffnessMatrix();
        end

        function [xF,errorVec] = compute(obj,fun,quadOrder)
            xF = copy(obj.trial);   
            obj.computeInitialGuess(fun,quadOrder);
            obj.updateDotProduct(obj.trial); 
            %obj.createRHSShapeDerivative(quadOrder);
            obj.chiN = obj.createRHSShapeFunction(fun,quadOrder);
            iter = 1;
            dJ0 = obj.computeCostGradient(quadOrder);
            error0 = Norm(dJ0,'L2');
            error  = inf;
            errorVec = [];
            while error >= 1e-6  && iter<=1000
                valOld = obj.trial.fValues;
                isAcceptable = false;
                while not(isAcceptable)
                    obj.createRHSDirectionalDerivative(quadOrder);
                    obj.solveProblem();
                    obj.updateDotProduct(obj.trial);
                    dJ = obj.computeCostGradient(quadOrder);
                    error = Norm(dJ,'L2');
                    if error <= 1.05*error0
                        obj.lineSearch = obj.lineSearch*1.1;
                        isAcceptable = true;
                    else
                        obj.lineSearch = obj.lineSearch/2;
                        obj.trial.setFValues(valOld);
                        obj.updateDotProduct(obj.trial);
                    end
                end
                iter = iter + 1;
                error0 = error;
                errorVec = [errorVec;error];
            end
            xF.setFValues(obj.trial.fValues);
        end

        function updateAlpha(obj,a)
            obj.alphaEps = a;
            obj.alpha    = a*obj.epsilon;
        end

        function updateBeta(obj,b)
            obj.betaEps = b;
            obj.beta    = b*obj.epsilon;
        end

        function updateEpsilon(obj,eps)
            obj.alpha      = obj.alphaEps*eps;
            obj.beta       = obj.betaEps*eps;
            obj.lineSearch = obj.lineSearch*(obj.epsilon/eps);
            obj.epsilon    = eps;
            obj.filter.updateEpsilon(max(obj.alpha,obj.beta));
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial      = LagrangianFunction.create(cParams.mesh, 1, 'P1');
            obj.mesh       = cParams.mesh;
            obj.epsilon    = obj.mesh.computeMeanCellSize();
            obj.alphaEps   = cParams.alpha;
            obj.betaEps    = cParams.beta;
            obj.alpha      = cParams.alpha*obj.epsilon;
            obj.beta       = cParams.beta*obj.epsilon;
            %obj.ub         = cParams.ub;
            %obj.lb         = cParams.lb;
            obj.lineSearch = 10;
        end

        function createFilterInitialGuess(obj)
            s.filterType = 'PDE';
            s.mesh  = obj.mesh;
            s.trial = LagrangianFunction.create(obj.mesh,1,obj.trial.order);
            f = Filter.create(s);
            f.updateEpsilon(max(obj.alpha,obj.beta));
            obj.filter = f;
        end

        function createDirection(obj,s)
            th = s.theta;            
            k  = [cosd(th);sind(th)];
            obj.direction = ConstantFunction.create(k,obj.mesh);
        end

        function createMassMatrix(obj)
            obj.M = IntegrateLHS(@(v,u) DP(v,u),obj.trial,obj.trial,obj.mesh,'Domain',2);
        end

        function createDirectionalStiffnessMatrix(obj)
            k       = obj.direction.constant;
          %  s.type  = 'AnisotropicStiffnessMatrix';
          %  s.mesh  = obj.mesh;
          %  s.test  = obj.trial;
          %  s.trial = obj.trial;
            A = ConstantFunction.create(k*k',obj.mesh);
          %  s.A     = A;
          %  LHS     = LHSIntegrator.create(s);
          %  K1   = LHS.compute();

            vF  = obj.trial;
            uF =  obj.trial;
            K2  = IntegrateLHS(@(u,v) DP(Grad(v),DP(A,Grad(u))'),vF,uF,obj.mesh,'Domain'); 
                
            obj.K = K2;

        end

        % function createRHSShapeDerivative(obj,quadOrder)
        %     s.mesh = obj.mesh;
        %     s.type = 'ShapeDerivative';
        %     s.quadratureOrder = quadOrder;
        %     s.test = obj.trial;
        %     obj.dNIntegrator = RHSIntegrator.create(s);
        % 
        % 
        % 
        % end

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
            s = obj.sVar;
            a = obj.alpha;
            b = obj.beta;
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            h = a.*maxFun - b.*minFun;
        end

        function createRHSDirectionalDerivative(obj,quadOrder)
            g         = obj.computeMeasureGradient();
            f         = obj.direction;
            %obj.kdhdN = obj.dNIntegrator.compute(f.*g);

            obj.kdhdN = IntegrateRHS(@(v) DP(Grad(v),f.*g),obj.trial,obj.mesh,quadOrder);

            
        end

        function g = computeMeasureGradient(obj)
            s      = obj.sVar;
            a      = obj.alpha;
            b      = obj.beta;
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            g      = s./(2*obj.lineSearch)-(a^2).*maxFun-(b^2).*minFun;
        end

        function solveProblem(obj)
            tau  = obj.lineSearch;
            LHS  = obj.M + (1/(2*tau))*obj.K;
            RHS  = obj.kdhdN + obj.chiN;
            x    = LHS\RHS;
            %rhoi = max(obj.lb,min(obj.ub,x));
            rhoi = x;
            obj.trial.setFValues(full(rhoi));
            %obj.isBoundFree = find((x-rhoi)==0);
        end

        function dJ = computeCostGradient(obj,quadOrder)
            dJUnc = obj.computeUnconstrainedCostGradient(quadOrder);
%             dj    = copy(dJUnc);
%             djV   = dj.fValues;
%             djV(obj.isBoundFree) = 0;
%             dj.setFValues(djV);
%             mu = obj.createDomainMax(-dj);
%             l  = obj.createDomainMax(dj);
%             dJ = dJUnc + mu - l; Future: box constraints
            dJ = dJUnc;
        end

        function dJ = computeUnconstrainedCostGradient(obj,quadOrder)
            [rhs1,rhs2,rhs3] = obj.computeCostGateaux(quadOrder);
            Ms   = obj.M;
            Ms   = diag(sum(Ms,1));
            fVal = Ms\(rhs1+rhs2+rhs3);
            dJ   = copy(obj.trial);
            dJ.setFValues(full(fVal));
        end

        function [rhs1,rhs2,rhs3] = computeCostGateaux(obj,quadOrder)
            a      = obj.alpha;
            b      = obj.beta;
            s      = obj.sVar;
            k      = obj.direction;
            maxFun = obj.createDomainMax(s);
            minFun = obj.createDomainMin(s);
            rhs1   = obj.createRHSShapeFunction(obj.trial,quadOrder);
            rhs2   = -obj.chiN;
            int3   = (a^2.*maxFun + b^2.*minFun).*k;
            rhs3   = IntegrateRHS(@(v) DP(Grad(v),int3),obj.trial,obj.mesh,quadOrder);
        end

        function RHS = createRHSShapeFunction(obj,fun,quadType)
            f   = @(v) DP(fun,v);
            RHS = IntegrateRHS(f,obj.trial,obj.trial.mesh,quadType);
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