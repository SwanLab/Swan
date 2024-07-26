classdef MicroComputer < handle

    properties (Access = public)
        computation
        variables
        solverType, solverMode
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = MicroComputer(cParams)
            obj.testName = cParams.testName;
            if isfield(cParams, 'solverType')
                obj.solverType = cParams.solverType;
            else
                obj.solverType = 'MONOLITHIC';
            end
            if isfield(cParams, 'solverMode')
                obj.solverMode = cParams.solverMode;
            else
                obj.solverMode = 'FLUC';
            end
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
            s.solverType = obj.solverType;
            s.solverMode = obj.solverMode;
            femSolver = PhysicalProblem.create(s);
            femSolver.solve();
            obj.computation = femSolver;
            obj.variables.Chomog = femSolver.Chomog;
        end

    end

end
