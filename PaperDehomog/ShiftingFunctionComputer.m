classdef ShiftingFunctionComputer < handle

    properties (Access = private)
        LHS
        RHS
    end
    
    properties (Access = private)
        mesh
        corrector
        interpolator
        test
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
            u = reshape(full(u),obj.mesh.nnodeElem,[]); 
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
            obj.test  = LagrangianFunction.create(obj.mesh,1,'P1D');

        end
         
        function computeLHS(obj)
            K     = IntegrateLHS(@(u,v) DP(Grad(u),Grad(v)),obj.test,obj.test,obj.mesh,2);
            In    = obj.interpolator;
            K = In'*K*In;
            obj.LHS = K;
        end

        function computeRHS(obj)
            gradC = Grad(obj.corrector);
            rhsF = IntegrateRHS(@(v) DP(Grad(v),gradC),obj.test,obj.mesh,2); 
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