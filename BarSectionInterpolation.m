classdef BarSectionInterpolation < handle
    
    properties (Access = public)
       sectionArea
       sectionInertia
       sectionAreaDerivative
       sectionInertiaDerivative
    end

    properties (Access = private)
        designVariable
        nVar
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
            var    = obj.designVariable.value;
            varNum = length(var)/2;
            r      = var(1:varNum);
            t      = var(varNum+1:end);
            A      = 4*pi.*r.*t;
            obj.sectionArea = A;
        end

        function computeSectionInertia(obj)
            var    = obj.designVariable.value;
            varNum = length(var)/2;
            r      = var(1:varNum);
            t      = var(varNum+1:end);
            I      = 8*pi*r*t*(r^2 + t^2);
            obj.sectionInertia = I;
        end

        function computeSectionAreaDerivative(obj)
            var = obj.designVariable.value;
            xR = var(1:obj.nVar);
            xT = var(obj.nVar+1:end);
            gR = 4*pi.*xT;
            gT = 4*pi.*xR;
            obj.sectionAreaDerivative = [gR;gT];   
        end

        function computeSectionInertiaDerivative(obj)
            xR = obj.designVariable(1:obj.nVar);
            xT = obj.designVariable(obj.nVar+1:end);
            gR = 8*pi.*xT*(3.*xR.^2 + xT.^2);
            gT = 8*pi.*xR*(3.*xT.^2 + xR.^2);
            obj.sectionInertiaDerivative = [gR;gT];
        end

    end
       
    methods (Access = private)
        
        function init(obj,designVar)
            obj.designVariable = designVar;
            obj.nVar = length(obj.designVariable.value)/2;
        end

    end
    
end