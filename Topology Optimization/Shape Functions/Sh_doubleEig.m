classdef Sh_doubleEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        eigNum
        eigModes
    end

    methods (Access = public)
        
        function obj = Sh_doubleEig(cParams)
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
            obj.eigNum = cParams.eigNum;
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