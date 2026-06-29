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
            traction = obj.tractionSeparation.computeFunction(obj.jump.fun);           
            F = IntegrateRHS(@(v) v'*Expand(traction,2),obj.jump,obj.cohesiveMesh.mesh,'Domain',quadOrder);
        end

        function [Ksec,Ktan] = computeDerivativeResidual(obj,u,quadOrder)
            obj.jump.updateJumpValues(u);
            Ksec = obj.computeDerivativeSecant(obj.jump.fun,quadOrder);
            Ktan = obj.computeDerivativeTangent(obj.jump.fun,quadOrder);
        end

        function updateDamageOld(obj,u)
            obj.jump.updateJumpValues(u);
            obj.tractionSeparation.updateDamageOld(obj.jump.fun);
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

        function LHS = computeLHS(obj,deriv,quadOrder)
            LHS = IntegrateLHS(@(u,v) DP(v', (DP(deriv, u'))'),obj.jump,obj.jump,obj.cohesiveMesh.mesh,'Domain',quadOrder);
        end

        function sec = computeDerivativeSecant(obj,jump,quadOrder)
            dt = obj.tractionSeparation.computeDerivativeSecant(jump);
            sec = obj.computeLHS(dt,quadOrder);
        end
        
        function tan = computeDerivativeTangent(obj,jump,quadOrder) 
            dt = obj.tractionSeparation.computeDerivativeTangent(jump);
            tan = obj.computeLHS(dt,quadOrder);
        end

    end
    
end