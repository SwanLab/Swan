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
            obj.a2  = obj.createConstantFunction(obj.alpha^2);
            obj.b2  = obj.createConstantFunction(obj.beta^2);
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);            
            obj.createRHSChi(fun,quadOrder);
            iter = 1;
            tolerance = 1;
            while tolerance >= 1e-4 
                obj.rhoOld.fValues = obj.trial.fValues;
                obj.createRHSDirectionalDerivative(quadOrder);
                obj.solveProblem();
                obj.updateDotProductPreviousGuess();
                obj.rhoDif.fValues = obj.rhoOld.fValues - obj.trial.fValues;
                tolerance = obj.rhoDif.computeL2norm()/obj.trial.computeL2norm();
                iter = iter + 1;
%                disp(iter);  
%                 disp(tolerance);
             end
           
%             obj.trial.plot
            xF.fValues = obj.trial.fValues;

        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial = LagrangianFunction.create(cParams.mesh, 1, 'P1'); % rho_eps
            obj.rhoOld = LagrangianFunction.create(cParams.mesh, 1, 'P1');
            obj.rhoDif = LagrangianFunction.create(cParams.mesh, 1, 'P1');
            obj.mesh  = cParams.mesh;
            obj.theta = cParams.theta;
            obj.alpha = cParams.alpha;
            obj.beta  = cParams.beta;
            obj.lineSearch = 100;
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
            LHS     = LHSintegrator.create(s);
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

        function createRHSDirectionalDerivative(obj,quadOrder)
            s.mesh = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSintegrator.create(s);
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
            m           = DomainFunction(s);
        end

        function m = createDomainMin(obj,sFun)
            s.operation = @(xV) min(zeros(size(xV(1,:,:))),sFun.evaluate(xV));
            m           = DomainFunction(s);
        end

        function solveProblem(obj)
            tau = obj.lineSearch;
            LHS = obj.M + (1/(2*tau))*obj.K;
            RHS = obj.rhsDer + obj.intChi;
            rhoi = LHS\RHS;
            obj.trial.fValues = rhoi;
        end

        function aF = createConstantFunction(obj,c)
            s.ndimf = length(c);
            s.operation = @(xV) obj.createOperation(c,xV);
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