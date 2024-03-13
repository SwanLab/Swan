classdef AcademicCost < handle
    
    properties (Access = private)
        costFunction
        gradientFunction
    end
    
    methods (Access = public)
        
        function obj = AcademicCost(cParams)
            obj.init(cParams)
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            J          = obj.costFunction(x.value);
            dJ.fValues = obj.gradientFunction(x.value);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.costFunction     = cParams.cF;
            obj.gradientFunction = cParams.gF;
        end
        
    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'AcademicCost';
        end
    end

end