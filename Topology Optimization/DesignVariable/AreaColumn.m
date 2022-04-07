classdef AreaColumn < DesignVariable
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        nElem
    end
    
    methods (Access = public)
        
        function obj = AreaColumn(cParams)
            obj.init(cParams)            
        end

        function A = getColumnArea(obj)
           x = obj.value;
           N = obj.nElem;
           A = x(1:N);
        end

        function V = computeVolum(obj)
            N = obj.getNelem();
            A = obj.getColumnArea();
            V = (1/N)*sum(A);
        end
        
        function gamma = getFirstEigenMode(obj)
           x = obj.value;
           N = obj.nElem;
           gamma = x(N+1);  
        end

        function nElem = getNelem(obj)
            nElem = obj.nElem;
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.value;
        end        
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
           obj.nElem = cParams.nElem; 
        end
        
    end
    
end