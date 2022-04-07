classdef Sh_doubleFirstEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        eigNum
        designVariable
    end
    
    properties (Access = private)
        eigModes
    end   
    
    methods (Access = public)
        
        function obj = Sh_doubleFirstEig(cParams)
            obj.init(cParams)
        end

        function computeFunctionAndGradient(obj,iter)
            obj.computeFunction(iter);
            obj.computeGradient();
        end
        
    end
    
    methods (Access = public)

        function computeFunction(obj,iter)
            eigN = obj.eigNum;
            fx = obj.eigModes.provideFunction(iter,eigN);
            obj.value = fx;
        end

        function computeGradient(obj,iter)
           eigN = obj.eigNum;
           dfdx = obj.eigModes.provideDerivative(iter,eigN);
           obj.gradient = dfdx';
        end

    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
            obj.eigModes = cParams.settings.eigMod;
            obj.eigNum = 1;
        end

  

    end    
    
end