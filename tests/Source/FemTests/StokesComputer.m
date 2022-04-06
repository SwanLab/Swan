classdef StokesComputer < handle

    properties (Access = public)
        computation
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = StokesComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            femSolver = StokesFEM.create(obj.testName);
            femSolver.computeVariables;
            obj.computation = femSolver;
        end
    end

end

