classdef ProblemSolver < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)

        function obj = ProblemSolver(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Static, Access = public)

        function solve(problem)
            [LHS, RHS, bc] = problem.getSolverData();
            Kred = bc.fullToReducedMatrix(LHS);
            Fred = bc.fullToReducedVector(RHS);
            u = Direct_Solver.solve(Kred,Fred);
            u = bc.reducedToFullVector(u);
            problem.setVariable(u);
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
        end

    end

end

