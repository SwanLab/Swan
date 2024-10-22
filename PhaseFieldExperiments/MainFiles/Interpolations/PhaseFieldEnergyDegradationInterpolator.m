classdef PhaseFieldEnergyDegradationInterpolator < handle

    properties (Access = public)
        fun
        dfun
        ddfun
    end

    properties (Access = private)
    end

    methods (Access = public)

        function obj = PhaseFieldEnergyDegradationInterpolator()
            obj.init()    
            obj.createEnergyDegradationFunctionAndDerivatives();
        end

    end

    methods (Access = private)

        function init(obj)  
        end

        function createEnergyDegradationFunctionAndDerivatives(obj)
            
            obj.fun      = @(phi) ((1-phi).^2);
            obj.dfun     = @(phi) -2.*(1-phi);
            obj.ddfun    = @(phi) 2;

            % obj.fun      = @(phi) ((1-sqrt(phi)).^2);
            % obj.dfun     = @(phi) (sqrt(phi)-1)./sqrt(phi);
            % obj.ddfun    = @(phi) 1/(2.*phi.^(3/2));

            % p = 4;
            % e = 1e-12;
            % obj.fun      = @(phi) ((1-phi.^p).^2);
            % obj.dfun     = @(phi) (-2*p).*(phi.^(p-1)).*(1-phi.^p);
            % obj.ddfun    = @(phi) (2*p^2).*((phi+e).^(2*(p-1))) + ((2*(p-1)*p).*((phi+e).^(p-2))).*((phi+e).^p-1);
        end

    end

end