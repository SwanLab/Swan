classdef Sh_doubleSecondEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        eigNum
    end

    properties (Access = private)
        eigModes
        bendingMat
        stiffnessMat
    end

    methods (Access = public)
        
        function obj = Sh_doubleSecondEig(cParams)
            obj.init(cParams)
        end

        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end
        
    end

    methods (Access = protected)

        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
            obj.eigModes = cParams.eigModes;
            obj.eigNum = 2;
        end
    end
    
    methods (Access = public)

        function computeFunction(obj)
           eigN = obj.eigNum;
           obj.value = obj.eigModes.provideFunction(eigN);
        end

        function computeGradient(obj)
            eigN = obj.eigNum;
            dfdx = obj.eigModes.provideDerivative(eigN);
            obj.gradient = dfdx';
        end
    end
    
end