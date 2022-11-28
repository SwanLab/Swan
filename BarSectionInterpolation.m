classdef BarSectionInterpolation < handle
    
    properties (Access = public)
       sectionArea
       sectionInertia
    end

    properties (Access = private)
        designVariable
    end

    methods (Access = public)
        
        function obj = BarSectionInterpolation(designVar)
            obj.init(designVar);
        end

        function computeSectionAreaAndInertia(obj)
            obj.computeSectionArea();
            obj.computeSectionInertia();
        end

        function computeSectionArea(obj)
            var    = obj.designVariable;
            varNum = length(var)/2;
            r      = var(1:varNum);
            t      = var(varNum+1:end);
            A      = 4*pi*r*t;
            obj.sectionArea = A;
        end

        function computeSectionInertia(obj)
            var    = obj.designVariable;
            varNum = length(var)/2;
            r      = var(1:varNum);
            t      = var(varNum+1:end);
            I      = 8*pi*r*t*(r^2 + t^2);
            obj.sectionInertia = I;
        end

    end
       
    methods (Access = private)
        
        function init(obj,designVar)
            obj.designVariable = designVar;
        end

    end
    
end