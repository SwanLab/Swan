classdef NonLinearFilterCircle < handle
    
    properties (Access = private)
        mesh
        trial
        q
        epsilon
    end

    properties (Access = private)
        M1
        M2
        RHSChi
        Kq
        RHS2
    end

    methods (Access = public)
        function obj = NonLinearFilterCircle(cParams)
            obj.init(cParams);
            obj.createMassMatrixFirstDirection();
            obj.createMassMatrixSecondDirection();
        end

        function xF = compute(obj,fun,quadOrder)
            xF = LagrangianFunction.create(obj.mesh, 1, obj.trial.order);
            obj.createRHSChiFirstDirection(fun,quadOrder);
            iter = 1;
            tolerance = 1;
            while tolerance >= 1e-5 && iter <= 5
                oldRho = obj.trial.fValues;
                obj.createKqFirstDirection(quadOrder);
                obj.solveFirstDirection();
                obj.createRHSSecondDirection(quadOrder);
                obj.solveSecondDirection();
                tolerance = norm(obj.trial.fValues - oldRho);
                iter = iter + 1;

                % Monitoring
%                 obj.trial.plot
%                 obj.q.plot
%                 gRho = Grad(obj.trial);
%                 energy_qDRho = gRho.*obj.q;
%                 E = Integrator.compute(energy_qDRho,obj.mesh,2);
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
            obj.trial   = LagrangianFunction.create(cParams.mesh, 1, 'P1'); % rho_eps
            obj.q       = LagrangianFunction.create(cParams.mesh, 2, 'P0'); % 2 = geom dim
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

        function createRHSChiFirstDirection(obj,fun,quadOrder)
            s.mesh     = obj.mesh;
            s.type     = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(fun,test);
            obj.RHSChi   = rhs;
        end


        function createKqFirstDirection(obj, quadOrder) % MISTAKE MAY BE HERE
            s.mesh = obj.mesh;
            s.type     = 'ShapeDerivative';
            s.quadratureOrder = quadOrder;
            int        = RHSintegrator.create(s);
            test       = obj.trial;
            rhs        = int.compute(obj.q,test);
            obj.Kq = rhs;
        end

        function createRHSSecondDirection(obj,quadOrder)
            s.mesh = obj.mesh;
            s.type = 'ShapeFunction';
            s.quadType = quadOrder;
            int        = RHSintegrator.create(s);
            nablaRho   = Grad(obj.trial);
            test       = obj.q;
            rhs        = int.compute(nablaRho,test);
            obj.RHS2   = -obj.epsilon^2*rhs;
        end



        function solveFirstDirection(obj)
            LHS = obj.M1;
            RHS = obj.RHSChi + obj.Kq;
            rhoi = LHS\RHS; % hard coded direct solver
            obj.trial.fValues = rhoi;
        end

        function solveSecondDirection(obj)
            LHS = obj.M2;
            RHS = obj.RHS2;
            qi = LHS \ RHS;
            obj.q.fValues = reshape(qi,[2,obj.mesh.nelem])'; % 2 = geo dim; P0 q obj.mesh.nelem / P1 q obj.mesh.nnodes
        end


        function itHas = hasEpsilonChanged(obj,eps)
            var   = abs(eps - obj.epsilon)/eps;
            itHas = var > 1e-15;
        end
    end
end