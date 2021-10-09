classdef StokesTest < PrecomputedVariableTest

    methods (Access = protected)
        function computeVariableThroughFemSolver(obj)
            femSolver = FEM.create(obj.testName);
            femSolver.computeVariables;
            obj.fem = femSolver;
        end        
    end

end

