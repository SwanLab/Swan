classdef QuadrilateralSection < SectionVariablesComputer
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = QuadrilateralSection(cParams)
            obj.init(cParams)
            
        end
        
        function I = computeInertia(obj)
            switch obj.designVariable.type
                case 'SquareColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    C = x(1:N);
                    I = (C.^4)/12;
                case 'RectangularColumn'
                    
                otherwise 
                    disp('Wrong Design Variable type');
            end
        end
        
        function dI = computeInertiaDerivative(obj)
            switch obj.designVariable.type
                case 'SquareColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    C = x(1:N);
                    dI = (4*C.^3)/12; 
                case 'RectangularColumn'
            end
        end
        
        function A = computeArea(obj)
            switch obj.designVariable.type
                case 'SquareColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    C = x(1:N);
                    A = C.^2;
                case 'RectangularColumn'
                    x = obj.designVariable.value;
                    N = mesh.nelem;
                    A = zeros(N,1);
                    for iElem=1:N
                        A(iElem)=x(iElem)*x(N+iElem);
                    end
            end
        end

        function dA = computeAreaDerivative(obj)
            switch obj.designVariable.type
                case 'SquareColumn'
                    x = obj.designVariable.value;
                    N = obj.mesh.nelem;
                    C = x(1:N);
                    dA = 2*C;
                case 'RectangularColumn'
            end
        end
        
    end
    
    methods (Access = private)
        
        
    end
    
end