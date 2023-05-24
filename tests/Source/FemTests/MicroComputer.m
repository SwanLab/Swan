classdef MicroComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = MicroComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            a.fileName = obj.testName;
            s = FemDataContainer(a);
%             femSolver = ElasticProblemMicro(s);
%             femSolver = NewElasticProblemMicro(s);
            femSolver = ElasticProblemMicro_Fast(s);
            femSolver.solve();
            obj.computation = femSolver;
            obj.variables.Chomog = femSolver.Chomog;
        end

    end

end
