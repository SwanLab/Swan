classdef NewFemSolver < NewTestComputationHandler

    methods (Static, Access = public)
        function solucio = compute(RHS, LHS)
            obj.fem = FEM.create(obj.testName);
            obj.createMaterialProperties();
            obj.fem.computeVariables();
        end
    end

end