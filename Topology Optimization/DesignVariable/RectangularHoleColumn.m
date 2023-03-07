classdef RectangularHoleColumn < DesignVariable
    
    properties (Access = public)
        nDesignVar = 4
        sectionType = 'Quadrilateral'
    end
    
    properties (Access = private)
        initValue
        initValueType
    end
       
    methods (Access = public)
        
        function obj = RectangularHoleColumn(cParams)
            obj.init(cParams);
            obj.type = 'RectangularHoleColumn';
            obj.createInitialValue();
        end
        
        function gamma = getFirstEigenMode(obj)
           x = obj.value;
           N = obj.mesh.nelem;
           nVar = obj.nDesignVar;
           gamma = x(nVar*N+1);  
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
                    b = ones(N,1);
                    c = 0.3*ones(N,1);
                    d = 0.3*ones(N,1);
                    x0 = [a;b;c;d;1];
                case 'Random'
                    a = rand(N,1);
                    b = rand(N,1);
                    c = 0.3*rand(N,1);
                    d = 0.3*rand(N,1);
                    x0 = [a;b;c;d;1];
                otherwise 
                    error('Invalid Initial Value Type.')
            end
            obj.update(x0);        
        end
    end

end