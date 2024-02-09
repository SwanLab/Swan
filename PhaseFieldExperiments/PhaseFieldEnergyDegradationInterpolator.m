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
            p = 2;

            obj.fun      = @(phi) ((1-phi).^p);
            obj.dfun     = @(phi) -p*((1-phi).^(p-1));
            obj.ddfun = @(phi) p*(p-1)*((1-phi).^(p-2));
        end

    end

end