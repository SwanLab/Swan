classdef Solver < handle
    properties
    end
    
    methods (Access = public, Static)
        function Sol = create(scale)
            switch scale
                case 'MACRO'                    
                    Sol = Solver_Dirichlet_Conditions;
                case 'MICRO'
                    Sol = Solver_Periodic;                    
            end
        end
    end
    
    methods (Abstract)
        solve(obj);
        setSolverVariables(obj);
    end
end

