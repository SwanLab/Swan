classdef ShiftingFunctionComputer < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
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
            s.fValues = u(:);
            s.order = 'P1D';
            sF = LagrangianFunction(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh     = cParams.mesh;
            obj.corrector    = cParams.corrector;
            obj.interpolator = cParams.interpolator;
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
            s.test  = LagrangianFunction.create(obj.mesh,1,'P1D');
            s.trial = LagrangianFunction.create(obj.mesh,1,'P1D');
            lhs = LHSIntegrator.create(s);
            K = lhs.compute();
        end

        function computeRHS(obj)
            gradC = Grad(obj.corrector);
            s.test = LagrangianFunction.create(obj.mesh,1,'P1D');
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            s.quadratureOrder = 2;
            rhs  = RHSIntegrator.create(s);
            rhsF = rhs.compute(gradC);
            In = obj.interpolator;
            rhsV = In'*rhsF;
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