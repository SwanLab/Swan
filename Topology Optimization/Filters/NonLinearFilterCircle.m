classdef NonLinearFilterCircle < handle
    
    properties (Access = private)
        mesh
        trial
        q
        epsilon
    end

    properties (Access = private)
        % Define new properties here
        M1
        M2
        RHS1
        Kq
        RHOi
        RHS2
    end

    methods (Access = public)
        function obj = NonLinearFilterCircle(cParams)
            obj.init(cParams);
            % Construct non-linear stuff...
            obj.createMassMatrixFirstDirection();
            obj.createMassMatrixSecondDirection();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            % solve non-linear filter...
            % let's start by creating a factory and define circle case (validation)
            obj.createRHSFirstDirection(fun,quadOrder);
            iter = 0;
            tolerance = 1;
            while tolerance >= 1e-5 && iter <= 100
                oldRho = obj.trial.fValues;
                obj.createKqFirstDirection(quadOrder);
                obj.solveFirstDirection();
                obj.createRhoiSecondDirection(quadOrder);
                obj.solveSecondDirection();
                tolerance = norm(obj.trial.fValues - oldRho);
                iter = iter + 1;

            end


            xF.fValues  = obj.trial.fValues;
        end

        function obj = updateEpsilon(obj,epsilon)
            if obj.hasEpsilonChanged(epsilon)
                obj.epsilon = epsilon;
                obj.computeLHS();
            end
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.trial   = LagrangianFunction.create(cParams.mesh, 1, cParams.trial.order); % rho_eps
            obj.q       = LagrangianFunction.create(cParams.mesh, 2, cParams.trial.order); % 2 = geom dim
            obj.mesh    = cParams.mesh;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function createMassMatrixFirstDirection(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.M1   = LHS.compute();
        end

        function createMassMatrixSecondDirection(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.q;
            s.trial = obj.q;
            s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.M2   = LHS.compute();
        end

        function createRHSFirstDirection(obj,fun,quadOrder)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(fun,test);
            obj.RHS1   = rhs;
        end


        function createKqFirstDirection(obj, quadOrder) % MISTAKE IS HERE
            s.mesh = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(obj.q,test); %Mass Matrix??
            obj.Kq = rhs;
        end

        function createRhoiSecondDirection(obj,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            nablaRho   = Grad(obj.trial);
            test       = obj.q;
            rhs        = int.compute(nablaRho,test);
            obj.RHS2   = -rhs;
        end



        function solveFirstDirection(obj)
            LHS = obj.M1;
            RHS = obj.RHS1 + obj.Kq;
            rhoi = LHS\RHS; % hard coded direct solver
            obj.trial.fValues = rhoi;
        end

        function solveSecondDirection(obj)
            LHS = obj.M2 .* (1/obj.epsilon^2);
            RHS = obj.RHS2;
            qi = LHS \ RHS; % hard coded direct solver
            obj.q.fValues = reshape(qi,[2,obj.mesh.nnodes])'; % 2 = geo dim
        end


        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end