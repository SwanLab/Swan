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
            [LHS, RHS, bc, newbc] = problem.getSolverData();
            LHSred = BCApplier.reduceLHS(LHS, newbc);
            RHSred = BCApplier.reduceLHS(RHS, newbc);
            u = Direct_Solver.solve(LHSred,RHSred);
            u = bc.reducedToFullVector(u);
            problem.setVariable(u);
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
        end

    end

end

