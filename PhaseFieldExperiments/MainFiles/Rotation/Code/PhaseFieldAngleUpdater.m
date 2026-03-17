classdef PhaseFieldAngleUpdater < handle
    
    properties (Access = private)
        monitor
    end

    methods (Access = public)

        function obj = PhaseFieldAngleUpdater(cParams)
            obj.init(cParams);
        end

        function [theta] = update(obj,u,phi)

        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.monitor    = cParams.monitor;
        end

    end

end