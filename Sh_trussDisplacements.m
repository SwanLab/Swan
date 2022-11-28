classdef Sh_trussDisplacements < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        physicalProblem
        barLength
        designVariable
    end
    
    methods (Access = public)
        
        function obj = Sh_trussDisplacements(cParams)
            obj.init(cParams);
        end
        
        
        function computeFunction(obj)
            p = obj.physicalProblem;
            p.solve();
            u = p.displacement;
            obj.value =
        end

        function computeGradient(obj)            
            xR = obj.designVariable(1:nVar);
            xT = obj.designVariable(nVar+1:end);
            gR = 4*pi*xT;
            gT = 4*pi*xR;
            obj.gradient = [gR;gT];         
        end
        
        
    end
    
    methods (Access = private)

        function obj = init(obj,cParams)
            obj.physicalProblem = cParams.phyProb;
            obj.maxDisplacement = cParams.maxDisplacement;
            obj.varN            = length(obj.designVariable)/2;
        end

    end
    
end