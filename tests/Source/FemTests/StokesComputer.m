classdef StokesComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = StokesComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = StokesDataContainer(a);
            femSolver = FEM.create(s);
            femSolver.computeVariables;
            obj.computation = femSolver;
            obj.variables = obj.computation.variables;
        end
    end

end

