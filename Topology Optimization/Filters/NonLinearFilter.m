classdef NonLinearFilter < handle
    
    properties (Access = private)
        mesh
        trial
        q
        epsilon
    end

    properties (Access = private)
        % Define new properties here
        M
        RHS1
        Kq
    end

    methods (Access = public)
        function obj = NonLinearFilter(cParams)
            obj.init(cParams);
            % Construct non-linear stuff...
            obj.createMassMatrix();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            % solve non-linear filter...
            % let's start by creating a factory and define circle case (validation)
            obj.createRHSFirstDirection(fun,quadOrder);
            obj.createKqFirstDirection(fun,quadOrder);
            obj.solveFirstDirection();

            %xF.fValues  = ;
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
            obj.q       = LagrangianFunction.create(cParams.mesh, 1, cParams.trial.order);
            obj.mesh    = cParams.mesh;
            obj.epsilon = cParams.mesh.computeMeanCellSize();
        end

        function createMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = obj.trial;
            s.trial = obj.trial;
            s.quadratureOrder = 2;
            LHS     = LHSintegrator.create(s);
            obj.M   = LHS.compute();
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


function createKqFirstDirection(obj,fun, quadOrder)
            s.mesh = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSintegrator.create(s);
            test = obj.q;
            rhs        = int.compute(fun,test);
            obj.Kq = rhs;
        end


        function solveFirstDirection(obj)
            LHS = obj.M;
            RHS = obj.RHS1 + obj.Kq;
            rhoi = LHS\RHS; % hard coded direct solver
            obj.trial.fValues = rhoi;
        end

        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end