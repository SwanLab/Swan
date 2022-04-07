classdef Sh_doubleSecondEig < ShapeFunctional
    
    properties (Access = private)
        nElem
        eigNum
        designVariable
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

        function computeFunctionAndGradient(obj,iter)
            obj.computeFunction(iter);
            obj.computeGradient();
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
            obj.eigModes = cParams.settings.eigMod;
            obj.eigNum = 2;
        end
    end
    
    methods (Access = public)

        function computeFunction(obj,iter)
           eigN = obj.eigNum;
           obj.value = obj.eigModes.provideFunction(iter,eigN);
        end

        function computeGradient(obj,iter)
            eigN = obj.eigNum;
            dfdx = obj.eigModes.provideDerivative(iter,eigN);
            obj.gradient = dfdx';
        end
    end
    
end