classdef Sh_volumeColumn < ShapeFunctional
    
    properties (Access = public)
        
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function obj = Sh_volumeColumn(cParams)
            obj.init(cParams)

        end
        function computeFunctionAndGradient(obj)
            obj.computeFunction();
            obj.computeGradient();
        end

    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
        end

    end

    methods (Access = public)

        function computeFunction(obj)
            V = obj.designVariable.computeVolum();
            fx = V-1;
            obj.value = fx;
        end

        function computeGradient(obj)
            nElem = obj.designVariable.mesh.nelem;
            dfdx = zeros(1,nElem+1);
            dfdx(1,1:nElem)=(1/nElem)*ones(1,nElem);
            obj.gradient = dfdx;
        end

    end

end