classdef ShiftingFunctionComputer < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        meshDisc
        mesh
        corrector
        interp
        field
        interpolator
    end
    
    methods (Access = public)
        
        function obj = ShiftingFunctionComputer(cParams)
            obj.init(cParams);
            obj.createField();
        end

        function sF = compute(obj)
            obj.computeLHS();
            obj.computeRHS();
            uC = obj.solveSystem();
            In = obj.interpolator;
            u  = In*uC; 
            u = reshape(u,obj.mesh.nnodeElem,[]); 
            s.mesh = obj.mesh;
            s.fValues(1,:,:) = u;
            sF = P1DiscontinuousFunction(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.corrector    = cParams.corrector;
            obj.interpolator = cParams.interpolator;
            obj.meshDisc     = obj.mesh.createDiscontinuousMesh();
        end

        function createField(obj)
            s.mesh               = obj.meshDisc;
            s.ndimf              = 1;
            s.interpolationOrder = obj.mesh.interpolation.order;
            obj.field = Field(s);
        end
         
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            In = obj.interpolator;
            K = In'*K*In;
            obj.LHS = K;
        end
        
        function K = computeStiffnessMatrix(obj)
            % Should be a P1DiscontinuousFunction instead!
            a.mesh = obj.meshDisc;
            a.fValues = zeros(obj.meshDisc.nnodes, 1);
            f = P1Function(a);

            a.mesh = obj.meshDisc;
            a.fValues = zeros(1, obj.meshDisc.nnodeElem, obj.mesh.nelem);
            fD = P1DiscontinuousFunction(a);

            s.mesh = obj.meshDisc;
            s.type = 'StiffnessMatrixFun';
            s.fun  = f;
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end

        function computeRHS(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('QUADRATIC');
            cG = obj.corrector.computeGradient(q);
            fGauss = permute(cG.fValues,[1 3 2]);
            s.fType     = 'Gauss';
            s.fGauss    = fGauss;
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