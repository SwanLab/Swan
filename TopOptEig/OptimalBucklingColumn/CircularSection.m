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
                    A = obj.computeArea();
                    I = (A.^2)';
                case 'RadiusColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    R = x(1:N);
                    I = pi^2*R.^4;
                otherwise 
                    disp('Wrong Design Variable type');
            end
        end
        
        function dI = computeInertiaDerivative(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    A = obj.computeArea();
                    dI = 2*A;
                case 'RadiusColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    R = x(1:N);
                    dI = 4*pi^2*R.^3;
            end
        end
        
        function A = computeArea(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    A = x(1:N,1);
                case 'RadiusColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    R = x(1:N);
                    A = R.^2.*pi;
            end
        end

        function dA = computeAreaDerivative(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    dA = 1;
                case 'RadiusColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    R = x(1:N);
                    dA = 2*pi*R;
            end
        end
        
    end
    
    methods (Access = private)
%         function I = computeAreaInertia(obj)
%             A = obj.computeArea();
%             I = (A.^2)';
%         end
        
%         function I = computeRadiusInertia(obj)
%             R = obj.designVariable.getColumnRadius();
%             I = pi^2*R^4;
%         end
        
%         function dI = computeAreaInertiaDerivative(obj)
%             A = obj.computeArea();
%             dI = 2*A;
%         end
%         
%         function dI = computeRadiusInertiaDerivative(obj)
%             R = obj.designVariable.getColumnRadius();
%             dI = 4*pi^2*R.^3;
%         end

%         function dA = computeRadiusAreaDerivative(obj)
%             R = obj.designVariable.getColumnRadius();
%             dA = 2*pi*R;
%         end
    end
    
end