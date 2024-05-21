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
            
            % obj.fun      = @(phi) ((1-phi.^2).^2);
            % obj.dfun     = @(phi) -4.*(phi-phi.^3);
            % obj.ddfun    = @(phi) -4 + 12.*phi.^2;
        end

    end

end