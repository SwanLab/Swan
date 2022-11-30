classdef Sh_trussTotalMass < handle
    
    properties (Access = public)
        value
        gradient
    end

    properties (Access = private)
        interpolator
        barLength
        designVariable
        varN
    end
    
    methods (Access = public)
        
        function obj = Sh_trussTotalMass(cParams)
            obj.init(cParams);
        end
        
        
        function computeFunction(obj)
            interp = obj.interpolator;
            interp.computeSectionArea();
            obj.value = sum(interp.sectionArea.*obj.barLength);
        end

        function computeGradient(obj)            
            interp = obj.interpolator;
            interp.computeSectionAreaDerivative();
            obj.gradient = interp.sectionAreaDerivative.*[obj.barLength;obj.barLength];
        end
        
        
    end
    
    methods (Access = private)

        function obj = init(obj,cParams)
            obj.interpolator   = cParams.interp;
            obj.designVariable = cParams.designVariable;
            obj.varN           = length(obj.designVariable)/2;
            obj.barLength      = cParams.barsLength;
        end

    end
    
end