classdef PhaseFieldFullQuadraticDissipation < handle

    properties (Access = public)
        interpolation
        constant
    end

    properties (Access = private)
        mesh
        xi
    end

    methods (Access = public)

        function obj = PhaseFieldFullQuadraticDissipation(cParams)
            obj.init(cParams)
            obj.computeDissipationFunctionAndDerivatives();
            obj.computeDissipationConstant()
        end
        
    end

    methods (Access = private)

        function init(obj,cParams)
            if cParams.xi < 0 || cParams.xi > 2
                error('Xi value incorrect [0,2]')
            else
                obj.xi = cParams.xi;
                obj.mesh = cParams.mesh;
            end
        end

        function computeDissipationConstant(obj)
            x = obj.xi;
            if x > 0 && x < 1
                obj.constant = (1/(1-x)^(3/2))*(0.5*(x^2)*log(x/(2*sqrt(1-x)+2-x)) +(2-x)*sqrt(1-x));
            elseif x == 1
                obj.constant = 8/3;
            else
                obj.constant = (1/(x-1)^(3/2))*(0.5*(x^2)*(pi/2 - asin((2-x)/x))-(2-x)*sqrt(x-1));
            end
        end

        function computeDissipationFunctionAndDerivatives(obj)
            s.operation = @(phi) obj.xi*phi + (1-obj.xi)*phi.^2;
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.interpolation.fun = DomainFunction(s);

            s.operation = @(phi) obj.xi + 2*(1-obj.xi)*phi;
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.interpolation.dfun = DomainFunction(s);

            s.operation = @(phi) 2*(1-obj.xi);
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.interpolation.ddfun = DomainFunction(s);

        end

    end

end