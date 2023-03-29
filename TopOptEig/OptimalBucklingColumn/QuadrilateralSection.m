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
                case 'RectangularHoleColumn'
                    [h,b,eh,eb]=obj.getFourValues();
                    I = (b.*h.^3-(b-2*eb).*(h-2*eh).^3)/12;
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
                case 'RectangularHoleColumn'
                    [h,b,eh,eb] = obj.getFourValues();
                    dIdh = (3*b.*h.^2-3*(b-2*eb).*(h-2*eh))/12;
                    dIdb = (h.^3-(h-2*eh).^3)/12;
                    dIdeh= ((b-2*eb).*(h-2*eh).^2)/2;
                    dIdeb= (h-2*eh).^3/2;
                    dI = [dIdh; dIdb; dIdeh; dIdeb];


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
                case 'RectangularHoleColumn'
                    [h,b,eh,eb] = obj.getFourValues();
                    h1 = h-2*eh;
                    b1 = b-2*eb;
                    A = (h.*b)-(h1.*b1);
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
                case 'RectangularHoleColumn'
                    [h,b,eh,eb] = obj.getFourValues();
                    h1 = h-2*eh;
                    b1 = b-2*eb;
                    dAdh = b-b1;
                    dAdb = h-h1;
                    dAdeh= 2*b1;
                    dAdeb= 2*h1;
                    dA = [dAdh; dAdb; dAdeh; dAdeb];
            end
        end
        
    end
    
    methods (Access = private)
        
        
    end
    
end