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
            % Should be a P1DiscontinuousFunction instead!
            a.mesh = obj.meshDisc;
            a.fValues = zeros(obj.meshDisc.nnodes, 1);
            f = P1Function(a);
            s.mesh = obj.meshDisc;
            s.type = 'StiffnessMatrixFun';
            s.fun  = f;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end

        function computeRHS(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            fG = obj.orientation.evaluate(q.posgp);            
            s.fType     = 'Gauss';
            s.fGauss    = fG;
            s.xGauss    = q.posgp;
            s.mesh      = obj.mesh;
            s.type      = obj.mesh.type;
            s.quadOrder = q.order;
            s.npnod     = obj.meshDisc.nnodes*1;
            s.type      = 'ShapeDerivative';
            s.globalConnec = obj.meshDisc.connec;
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute();
            In = obj.interpolator;
            rhsV = In'*rhsV;
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