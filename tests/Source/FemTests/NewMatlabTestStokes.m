classdef NewMatlabTestStokes < NewMatlabTest

    methods (Access = protected)
        
%         function selectComputedVar(obj)
%             obj.computedVar{1} = obj.fem.variables.u;
%             obj.computedVar{2} = obj.fem.variables.p;
%         end

    end

    methods (Access = protected) % heredat de testFemComputation
        function computeVariableThroughFemSolver(obj)
            femSolver = FEM.create(obj.testName);
            femSolver.computeVariables;
            obj.fem = femSolver;
        end        
    end

end

