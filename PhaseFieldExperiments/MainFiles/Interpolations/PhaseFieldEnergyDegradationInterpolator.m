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
            s.operation = @(phi) ((1-phi).^2);
            s.ndimf = 1;
            obj.fun = DomainFunction(s);

            s.operation = @(phi) -2.*(1-phi);
            s.ndimf = 1;
            obj.dfun = DomainFunction(s);

            s.operation = @(phi) 2;
            s.ndimf = 1;
            obj.ddfun = DomainFunction(s);
        end

    end

end