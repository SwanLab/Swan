classdef MicroComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = MicroComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            femSolver = Elastic_Problem_Micro.create(obj.testName);
            props.kappa = .9107;
            props.mu    = .3446;
            femSolver.setMatProps(props);
            femSolver.computeChomog;
            obj.computation = femSolver;
        end
    end

end

