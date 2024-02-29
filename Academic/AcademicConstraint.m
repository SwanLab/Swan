classdef AcademicConstraint < handle
    
    properties (Access = private)
        gradientFunction
        constraintFunction
    end
    
    methods (Access = public)
        
        function obj = AcademicConstraint(cParams)
            obj.init(cParams)
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            J          = obj.constraintFunction(x.value);
            dJ.fValues = obj.gradientFunction(x.value);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.constraintFunction = cParams.cF;
            obj.gradientFunction   = cParams.gF;
        end

    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'AcademicConstraint';
        end
    end

end