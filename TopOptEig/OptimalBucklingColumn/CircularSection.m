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
                    R = obj.getSingleValue();
                    I = pi^2*R.^4;
                case 'HoleColumn'
                    [r,e] = obj.getDoubleValue();
                    I = pi*((r+e).^4-r.^4)/4;
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
                    R = obj.getSingleValue();
                    dI = 4*pi^2*R.^3;
                case 'HoleColumn'
                    [r,e] = obj.getDoubleValue();
                    dIdr = pi*((r+e).^3-r.^3);
                    dIde = pi*(r+e).^3;
                    dI = [dIdr; dIde];
            end
        end
        
        function A = computeArea(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    A = obj.val;
                case 'RadiusColumn'
                    R = obj.getSingleValue();
                    A = R.^2*pi;
                case 'HoleColumn'
                    [r,e] = obj.getDoubleValue();
                    Rint = r;
                    Rext = r+e;
                    A = pi*(Rext.^2-Rint.^2);

            end
        end

        function dA = computeAreaDerivative(obj)
            switch obj.designVariable.type
                case 'AreaColumn'
                    dA = 1;
                case 'RadiusColumn'
                    R = obj.getSingleValue();
                    dA = 2*pi*R;
                case 'HoleColumn'
                    [r,e] = obj.getDoubleValue();
                    dAdr = 2*pi*e;
                    dAde = 2*pi*(r+e);
                    dA = [dAdr; dAde];
            end
        end
        
    end
    
    methods (Access = private)
        
    end
    
end