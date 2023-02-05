classdef SquareColumn < DesignVariable
    
    properties (Access = public)
        nDesignVar = 1
    end
    
    properties (Access = private)
        initValue
        initValueType
    end
    
    methods (Access = public)
        
        function obj = SquareColumn(cParams)
            obj.init(cParams)
            obj.type = 'SquareColumn';
            obj.createInitialValue();
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
                    x0 = sqrt(rand(N+1,1));
                case 'External Value'
                    x0 = obj.initValue;
                    x0=x0+norm(x0)*rand(size(x0))*0.01;
                otherwise 
                    error('Invalid Initial Value Type.')
            end
            obj.update(x0);        
        end
    end
    
end