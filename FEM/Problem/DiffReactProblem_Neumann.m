classdef DiffReactProblem_Neumann < DiffReactProblem

    methods (Access = public)
        
        function computeVariables(obj,x)
            bc   = obj.boundaryConditions;
            xRed = bc.fullToReducedVector(x);
            LHS  = obj.computeLHS();
            xReg = obj.solver.solve(LHS,xRed);
            obj.variables.x = bc.reducedToFullVector(xReg);
        end
        
        function LHS = computeLHS(obj)
            LHS = obj.epsilon^2*obj.K + obj.M;
            LHS = obj.boundaryConditions.fullToReducedMatrix(LHS);
        end
        
    end
    
end
