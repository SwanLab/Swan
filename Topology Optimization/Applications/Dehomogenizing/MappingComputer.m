classdef MappingComputer < handle

    properties (Access = private)
        LHS
        RHS
    end

    properties (Access = private)
        meshDisc
        mesh
        orientation
        interp
        interpolator
    end

    methods (Access = public)

        function obj = MappingComputer(cParams)
            obj.init(cParams);
        end

        function uF = compute(obj)
            obj.computeLHS();
            obj.computeRHS();
            uC = obj.solveSystem();
            In = obj.interpolator;
            u  = In*uC;
            uV(1,:,:) = reshape(u,obj.mesh.nnodeElem,[]);
            s.mesh    = obj.mesh;
            s.fValues = uV;
            uF = P1DiscontinuousFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation  = cParams.orientation;
            obj.interpolator = cParams.interpolator;
            obj.meshDisc     = obj.mesh.createDiscontinuousMesh();
        end
        
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            In = obj.interpolator;
            Kn = In'*K*In;
            obj.LHS = Kn;
        end

        function K = computeStiffnessMatrix(obj)
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            s.test  = P1DiscontinuousFunction.create(obj.mesh, 1);
            s.trial = P1DiscontinuousFunction.create(obj.mesh, 1);
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end

        function computeRHS(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
%             fG = obj.orientation.evaluate(q.posgp);
            s.mesh  = obj.meshDisc;
            s.type = 'ShapeDerivative';
            s.quadratureOrder = q.order;
            rhs  = RHSintegrator.create(s);
            rhsF = rhs.compute(obj.orientation,obj.orientation);
            In = obj.interpolator;
            rhsV = In'*rhsF.fValues;
            obj.RHS = rhsV;
        end

        function u = solveSystem(obj)
            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(obj.LHS,obj.RHS);
            u = u(1:end);
        end

    end

end