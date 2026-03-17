classdef PhaseFieldAngleUpdater < handle
    
    properties (Access = private)
    end

    methods (Access = public)

        function obj = PhaseFieldAngleUpdater()
            obj.init();
        end

        function [theta] = update(obj,u,theta,phi)
            theta.setFValues((0)*ones(size(theta.fValues)));
            %theta.setFValues((pi/2)*ones(size(theta.fValues)));
        end

    end

    methods (Access = private)

        function init(obj)
        end

    end

end