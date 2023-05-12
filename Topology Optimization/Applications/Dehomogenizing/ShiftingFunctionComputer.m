classdef ShiftingFunctionComputer < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        meshDisc
        mesh
        corrector
        interpolator
    end
    
    methods (Access = public)
        
        function obj = ShiftingFunctionComputer(cParams)
            obj.init(cParams);
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
         
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            In = obj.interpolator;
            K = In'*K*In;
            obj.LHS = K;
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
            cG = obj.corrector.computeGradient(q);

            s.mesh = obj.meshDisc;
            s.type = 'ShapeDerivative';
            s.quadratureOrder = 'QUADRATIC';
            rhs  = RHSintegrator.create(s);
            rhsF = rhs.compute(cG);
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