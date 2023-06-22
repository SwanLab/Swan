classdef HoleColumn < DesignVariable
    
    properties (Access = public)
        nDesignVar = 2
        sectionType = 'Circular'
    end
    
    properties (Access = private)
        initValue
        initValueType
    end
       
    methods (Access = public)
        
        function obj = HoleColumn(cParams)
            obj.init(cParams);
            obj.type = 'HoleColumn';
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
                    r = ones(N,1);
                    e = 0.3*ones(N,1);
                    x0 = [r;e;1];
                case 'Random'
                    r = rand(N,1);
                    e = 0.3*rand(N,1);
                    x0 = [r;e;1];
                case 'Sinus'
                    x = linspace(0,5*pi,N)';
                    r = 0.4*sin(x)+1;
                    e = 0.2*sin(x)+0.3;
                    x0 = [r;e;1];
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