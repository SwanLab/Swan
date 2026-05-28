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
            obj.jump.updateJumpValues(u);
            traction  = obj.tractionSeparation.computeFunction(obj.jump.fun);
            fun       = DP(traction, obj.jump.fun);
            E         = Integrator.compute(fun,obj.cohesiveMesh.mesh,quadOrder);
        end

        function F = computeResidual(obj,u,quadOrder)
            obj.jump.updateJumpValues(u);
            traction = obj.tractionSeparation.computeFunction(obj.jump.fun);           
            F = IntegrateRHS(@(v) v'*Expand(traction,2),obj.jump,obj.cohesiveMesh.mesh,'Domain',quadOrder);
        end

        function H = computeDerivativeResidual(obj,u,quadOrder)
            obj.jump.updateJumpValues(u);
            deriv = obj.tractionSeparation.computeDerivative(obj.jump.fun);
            H = IntegrateLHS(@(u,v) DP(v', (DP(deriv, u'))'),...
                obj.jump,obj.jump,obj.cohesiveMesh.mesh,'Domain',quadOrder);
        end

        function updateLambdaOld(obj)
            obj.tractionSeparation.updateLambdaOld();
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.cohesiveMesh     = cParams.cohesiveMesh;
            obj.material         = cParams.material;
            obj.tractionSeparation = cParams.tractionSeparation;
        end

        function createJumpFunction(obj)
            s.cohesiveMesh = obj.cohesiveMesh;
            s.ndimf = obj.cohesiveMesh.fullMesh.ndim;
            s.uFun  = LagrangianFunction.create(obj.cohesiveMesh.fullMesh,s.ndimf,'P1');
            obj.jump = Jump(s);
        end

    end
    
end