classdef PhaseFieldEnergyDegradationInterpolator < handle

    properties (Access = public)
        fun
        dfun
        ddfun
    end

    properties (Access = private)
        mesh
    end

    methods (Access = public)

        function obj = PhaseFieldEnergyDegradationInterpolator(cParams)
            obj.init(cParams)    
            obj.createEnergyDegradationFunctionAndDerivatives();
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end

        function createEnergyDegradationFunctionAndDerivatives(obj)
            s.operation = @(phi) ((1-phi).^2);
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.fun = DomainFunction(s);

            s.operation = @(phi) -2.*(1-phi);
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.dfun = DomainFunction(s);

            s.operation = @(phi) 2;
            s.ndimf = 1;
            s.mesh = obj.mesh;
            obj.ddfun = DomainFunction(s);
        end

    end

end