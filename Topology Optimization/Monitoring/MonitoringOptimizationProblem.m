classdef MonitoringOptimizationProblem < handle  % will be topologymonitoring

    properties (Access = private)
        cost
        functionals
        constraints
        dualVariables
        designVariable

    end

    methods
        function obj = MonitoringOptimizationProblem(inputArg1,inputArg2)
            %UNTITLED13 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end