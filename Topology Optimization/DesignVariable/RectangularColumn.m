classdef RectangularColumn < DesignVariable
    
    properties (Access = public)
        nDesignVar = 2
        setionType = 'Quadrilateral'
    end
    
    properties (Access = private)
        initValue
        initValueType
    end
       
    methods (Access = public)
        
        function obj = RectangularColumn(cParams)
            obj.init(cParams);
            obj.type = 'RectangularColumn';
            obj.createInitialValue();
        end
        
        function gamma = getFirstEigenMode(obj)
           x = obj.value;
           N = obj.mesh.nelem;
           gamma = x(2*N+1);  
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
                    a = ones(N,1);
                    b = 0.7*ones(N,1);
                    x0 = [a;b;1];
                case 'Random'
                    a = rand(N,1);
                    b = 0.7*rand(N,1);
                    x0 = [a;b;1];
                otherwise 
                    error('Invalid Initial Value Type.')
            end
            obj.update(x0);        
        end
    end

end