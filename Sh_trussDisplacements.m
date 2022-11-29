classdef Sh_trussDisplacements < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        physicalProblem
        designVariable
        normVal
        qVal
        maxDisplacement
    end
    
    methods (Access = public)
        
        function obj = Sh_trussDisplacements(cParams)
            obj.init(cParams);
        end
        
        
        function computeFunction(obj)
            p = obj.physicalProblem;
            p.solve();
            u = p.displacement;
            obj.value = norm(u,obj.qVal);
        end

        function computeGradient(obj)            
            p = obj.physicalProblem;
            p.solveDisplacementAdjoint();
            obj.gradient = -p.Adjoint.*p.stifnessDerivative;
        end
        
        
    end
    
    methods (Access = private)

        function obj = init(obj,cParams)
            obj.physicalProblem = cParams.phyProb;
            obj.maxDisplacement = cParams.maxDisplacement;
            obj.normVal         = cParams.normVal;
            obj.qVal            = cParams.qVal;
            obj.varN            = length(obj.designVariable)/2;
        end

    end
    
end