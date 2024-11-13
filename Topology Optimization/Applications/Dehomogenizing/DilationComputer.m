classdef DilationComputer < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        mesh
        orientationVector
        dilation
    end
    
    methods (Access = public)
        
        function obj = DilationComputer(cParams)
            obj.init(cParams);
            obj.createDilationFun();
        end

        function rF = compute(obj)
            obj.computeLHS();
            obj.computeRHS();
            r = obj.solveSystem();
            s.mesh = obj.mesh;
            s.fValues = r;
            s.order = 'P1';
            rF = LagrangianFunction(s);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh               = cParams.mesh;
            obj.orientationVector  = cParams.orientationVector;
        end
       
        function computeLHS(obj)
            K = obj.computeStiffnessMatrix();
            I = ones(size(K,1),1);
            obj.LHS = [K,I;I',0];
        end
        
        function K = computeStiffnessMatrix(obj)
            s.test  = obj.dilation;
            s.trial = obj.dilation;
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
        end

        function createDilationFun(obj)
            obj.dilation = LagrangianFunction.create(obj.mesh, 1, 'P1');
        end
        
        function computeRHS(obj)
            s.mesh      = obj.mesh;
            s.operation = @(x) obj.computeFieldTimesDivField(x);
            gradT       = DomainFunction(s);
            s.mesh = obj.mesh;
            s.type = 'ShapeDerivative';
            s.quadratureOrder = 3;
            test = LagrangianFunction.create(obj.mesh,1,'P1');
            rhs  = RHSintegrator.create(s);
            rhsV = rhs.compute(gradT,test);
            obj.RHS = [rhsV;0];
        end
        
        function gradT = computeFieldTimesDivField(obj,xV)
            a1    = obj.orientationVector{1};
            a2    = obj.orientationVector{2};
            aDa1  = a1.computeFieldTimesDivergence(xV);
            aDa2  = a2.computeFieldTimesDivergence(xV);  
            gradT = -aDa1.fValues - aDa2.fValues;
        end
        
        function u = solveSystem(obj)
            a.type = 'DIRECT';
            s = Solver.create(a);
            u = s.solve(obj.LHS,obj.RHS);
            u = u(1:end-1);
        end
        
    end
    
end