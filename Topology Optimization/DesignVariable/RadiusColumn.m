classdef RadiusColumn < DesignVariable
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        initValue
        initValueType
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = RadiusColumn(cParams)
            obj.init(cParams);
            obj.type = 'RadiusColumn';
            obj.createInitialValue();
        end
        
%         function R = getColumnRadius(obj)
%             x = obj.value;
%             N = obj.mesh.nelem;
%             R = x(1:N);
%         end
        
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
                    x0(1:N+1,1) = sqrt(1/pi);
                case 'Random'
                    x0 = sqrt(rand(N+1,1)/pi);
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