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
            obj.createJumpFunction(cParams);
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

            % func = Bc.' * deriv * Bc
            func = DP(Bc, DP(deriv,Bc));

            H = IntegrateLHS(@(u,v)  func...
                ,obj.test,obj.test,obj.cohesiveMesh.mesh,'Domain',quadOrder);
        
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.cohesiveMesh     = cParams.cohesiveMesh;
            obj.material = cParams.material;
            obj.test     = cParams.test;
            obj.tractionSeparation = cParams.tractionSeparation;
        end

        function createJumpFunction(obj,cParams)
            s.cohesiveMesh = cParams.cohesiveMesh;
            s.ndimf = cParams.ndimf;
            s.uFun = ConstantFunction.create(zeros(cParams.ndimf),cParams.cohesiveMesh.mesh);
            obj.jump = Jump(s);
        end
    end
    
end