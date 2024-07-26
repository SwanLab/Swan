classdef AcademicConstraint < handle
    
    properties (Access = private)
        constraintFunction
        gradientFunction
    end
    
    methods (Access = public)
        
        function obj = AcademicConstraint(cParams)
            obj.init(cParams)
        end
        
        function [J,dJ] = computeFunctionAndGradient(obj,x)
            J          = obj.constraintFunction(x.fun.fValues);
            dJ.fValues = obj.gradientFunction(x.fun.fValues);
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