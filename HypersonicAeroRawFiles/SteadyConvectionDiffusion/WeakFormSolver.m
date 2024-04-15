classdef WeakFormSolver < handle

    methods (Static, Access = public)
        function obj = create(cParams)
            w = WeakFormSolverFactory();
            obj = w.create(cParams);
        end
    end
end