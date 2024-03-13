classdef AreaColumn < DesignVariable
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
       
    methods (Access = public)
        
        function obj = AreaColumn(cParams)
            obj.init(cParams);
            obj.type = 'AreaColumn';
            obj.createInitialValue();
        end

        function A = getColumnArea(obj)
           x = obj.value;
           N = obj.mesh.nelem;
           A = x(1:N,1);
        end

        function V = computeVolum(obj)
            N = obj.mesh.nelem;
            A = obj.getColumnArea();
            V = (1/N)*sum(A);
        end
        
        function gamma = getFirstEigenMode(obj)
           x = obj.value;
           N = obj.mesh.nelem;
           gamma = x(N+1);
        end

        function v = getVariablesToPlot(obj)
            v{1} = obj.value;
        end

        function norm = computeL2normIncrement(obj)
            norm = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
        function createInitialValue(obj)
            N = obj.mesh.nelem;
            x0 = 1*rand(N+1,1);%ones(N+1,1);               
            obj.update(x0);        
        end
    end

end