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
                    C = obj.getSingleValue();
                    I = (C.^4)/12;
                case 'RectangularColumn'
                    [a,b]=obj.getDoubleValue();
                    I = a.*(b.^3)/12;
                otherwise 
                    disp('Wrong Design Variable type');
            end
        end
        
        function dI = computeInertiaDerivative(obj)
            switch obj.designVariable.type
                case 'SquareColumn'
                    C = obj.getSingleValue();
                    dI = (4*C.^3)/12; 
                case 'RectangularColumn'
                    [a,b]=obj.getDoubleValue();
                    dIda = b.^3/12;
                    dIdb = 3*a.*(b.^2)/12;
                    dI = [dIda; dIdb];
            end
        end
        
        function A = computeArea(obj)
            switch obj.designVariable.type
                case 'SquareColumn'
                    C = obj.getSingleValue();
                    A = C.^2;
                case 'RectangularColumn'
                    [a,b]=obj.getDoubleValue();
                    A = a.*b;
            end
        end

        function dA = computeAreaDerivative(obj)
            switch obj.designVariable.type
                case 'SquareColumn'
                    C = obj.getSingleValue();
                    dA = 2*C;
                case 'RectangularColumn'
                    [a,b]=obj.getDoubleValue();
                    dAda = b;
                    dAdb = a;
                    dA = [dAda; dAdb];
            end
        end
        
    end
    
    methods (Access = private)
        
        
    end
    
end