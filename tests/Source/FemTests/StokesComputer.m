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
            femSolver = PhysicalProblem.create(s);
            femSolver.computeVariables;
            obj.computation = femSolver;
            u = obj.computation.velocityFun.fValues;
            obj.variables.p = obj.computation.pressureFun.fValues;
            obj.variables.u = reshape(u', [numel(u) 1]);
        end
    end

end

