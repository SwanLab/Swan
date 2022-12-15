classdef CircularSection < SectionVariablesComputer
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = CircularSection(cParams)
            obj.init(cParams)
            
        end
        
        function I = computeInertia(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    I = obj.computeAreaInertia();
                case 'RadiusColumn'
                    I = obj.computeRadiusInertia();
                otherwise 
                    disp('Wrong Design Variable type');
            end
        end
        
        function dI = computeInertiaDerivative(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    dI = obj.computeAreaInertiaDerivative();
                case 'RadiusColumn'
                    dI = obj.computeRadiusInertiaDerivative();
            end
        end

        function dA = computeAreaDerivative(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    dA = 1;
                case 'RadiusColumn'
                    dA = obj.computeRadiusAreaDerivative();
            end
        end
        
    end
    
    methods (Access = private)
        function I = computeAreaInertia(obj)
            A = obj.designVariable.getColumnArea();
            I = (A.^2)';
        end
        
        function I = computeRadiusInertia(obj)
            R = obj.designVariable.getColumnRadius();
            I = pi^2*R^4;
        end
        
        function dI = computeAreaInertiaDerivative(obj)
            A = obj.designVariable.getColumnArea();
            dI = 2*A;
        end
        
        function dI = computeRadiusInertiaDerivative(obj)
            R = obj.designVariable.getColumnRadius();
            dI = 4*pi^2*R.^3;
        end

        function dA = computeRadiusAreaDerivative(obj)
            R = obj.designVariable.getColumnRadius();
            dA = 2*pi*R;
        end
    end
    
end