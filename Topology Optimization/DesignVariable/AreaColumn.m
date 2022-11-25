classdef AreaColumn < DesignVariable
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        initValue
        initValueType
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
            obj.initValue = cParams.initValue;
            obj.initValueType = cParams.initValueType;
        end
        
        function createInitialValue(obj)
            N = obj.mesh.nelem;
            switch obj.initValueType
                case 'Constant'
                    x0 = ones(N+1,1);
                case 'Random'
                    x0 = rand(N+1,1);
                case 'External Value'
                    x0 = obj.initValue;
                    x0=x0+norm(x0)*rand(size(x0))*0.01;
                otherwise 
                    error('Invalid Initial Value Type.')
            end
            %x0 = 2*rand(N+1,1);
            obj.update(x0);        
        end
    end

end