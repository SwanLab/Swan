classdef Sh_doubleFirstEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        eigNum
    end
    
    properties (Access = private)
        eigModes
    end   
    
    methods (Access = public)
        
        function obj = Sh_doubleFirstEig(cParams)
            obj.init(cParams)
        end

        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end
        
    end
    
    methods (Access = public)

        function computeFunction(obj)
            eigN = obj.eigNum;
            fx = obj.eigModes.provideFunction(eigN);
            obj.value = fx;
        end

        function computeGradient(obj)
           eigN = obj.eigNum;
           dfdx = obj.eigModes.provideDerivative(eigN);
           obj.gradient = dfdx';
        end

    end

    methods (Access = protected)
        
        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
            obj.eigModes = cParams.eigModes;
            obj.eigNum = 1;
        end

  

    end    
    
end