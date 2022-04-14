classdef DiffReactProblem_Neumann < DiffReactProblem

    methods (Access = public)
        
        function computeVariables(obj,x)
            bc   = obj.boundaryConditions;
            xRed = bc.fullToReducedVector(x);
            LHS  = obj.computeLHS();
            x = obj.solver.solve(LHS,xRed);
            obj.variables.x = bc.reducedToFullVector(x);
        end
        
        function LHS = computeLHS(obj)
            LHS = obj.epsilon^2*obj.K + obj.M;
            LHS = obj.boundaryConditions.fullToReducedMatrix(LHS);
        end
        
    end
    
end
