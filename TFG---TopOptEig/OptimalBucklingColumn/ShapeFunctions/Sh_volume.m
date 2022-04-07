classdef Sh_volume < ShapeFunctional
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        nElem
        designVariable
    end
    
    methods (Access = public)
        
        function obj = Sh_volume(cParams)
            obj.init(cParams)

        end
        function computeFunctionAndGradient(obj,iter)
            obj.computeFunction(iter);
            obj.computeGradient(iter);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nElem = cParams.nElem;
            obj.designVariable = cParams.designVariable;
        end

    end

    methods (Access = public)

        function computeFunction(obj,iter)
            x = obj.designVariable.value;
            N = obj.nElem;
            fx = (1/N)*sum(x(1:N))-1;
            obj.value = fx;
        end

        function computeGradient(obj,iter)
            N = obj.nElem;
            dfdx = zeros(1,N+1);
            dfdx(1,1:N)=(1/N)*ones(1,N);
            obj.gradient = dfdx;
        end

    end

end