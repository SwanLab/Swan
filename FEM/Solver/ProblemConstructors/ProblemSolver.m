classdef ProblemSolver < handle

    properties (Access = protected)
        boundaryConditions
        solver
    end

    methods (Access = public, Static)

        function obj = create(s)
            p = ProblemConstructorFactory();
            obj = p.create(s);
        end

    end
    
    methods (Access = public)

        function obj = ProblemSolver(cParams)
            obj.init(cParams);
        end

    end

    methods (Access = protected)

        function init(obj, cParams)
            obj.boundaryConditions = cParams.boundaryConditions;
            obj.solver             = cParams.solver;
        end

    end
    
    methods (Abstract)
        [u,L] = solve(obj,cParams)
    end

end