classdef CohesiveTermFunctional < handle
    
    properties (Access = public)
    end
    
    properties (Access = private)
        cohesiveMesh
        test
    end
    
    properties (Access = private)
        material
        jump
        tractionSeparation
    end
    
    methods (Access = public)
        
        function obj = CohesiveTermFunctional(cParams)
            obj.init(cParams)
            obj.createJumpFunction();
        end

        function E = computeCost(obj,u,quadOrder)

        end

        function F = computeResidual(obj,u,quadOrder)
            %F =  grau llibertat desplaçament 
        end

        function H = computeDerivativeResidual(obj,u,quadOrder)
            obj.jump.updateJumpValues(u);
            jumpFun = obj.jump.fun;
            deriv = obj.tractionSeparation.computeDerivative(jumpFun);
            func = DP(Bc, DP(deriv,Bc));            % func = Bc.' * deriv * Bc
            H = IntegrateLHS(@(u,v)  func...
                ,obj.test,obj.test,obj.cohesiveMesh.mesh,'Domain',quadOrder);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.cohesiveMesh     = cParams.cohesiveMesh;
            obj.material = cParams.material;
            obj.test     = cParams.cohTest;
            obj.tractionSeparation = cParams.tractionSeparation;
        end

        function createJumpFunction(obj)
            s.cohesiveMesh = obj.cohesiveMesh;
            s.ndimf = obj.cohesiveMesh.mesh.ndim;
            s.uFun  = LagrangianFunction.create(obj.cohesiveMesh.mesh,s.ndimf,'P1');
            obj.jump = Jump(s);
        end
    end
    
end