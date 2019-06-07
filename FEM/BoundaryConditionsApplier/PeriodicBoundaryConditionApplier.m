classdef PeriodicBoundaryConditionApplier < BoundaryConditionsApplier
    
    properties (Access = private)
        nfields
        dof        
    end
    
    methods (Access = public)
        
        function obj = PeriodicBoundaryConditionApplier(cParams)
            obj.nfields = cParams.nfields;
            obj.dof = cParams.dof;            
        end
        
        function Ared = fullToReducedMatrix(obj,A)                
            vF = obj.dof.free;
            vP = obj.dof.periodic_free;
            vQ = obj.dof.periodic_constrained;
            vI = setdiff(vF{1},vP);
            
            A_II = A(vI,vI);
            A_IP = A(vI,vP) + A(vI,vQ); %Grouping P and Q nodal values
            A_PI = A(vP,vI) + A(vQ,vI); % Adding P  and Q equation
            A_PP = A(vP,vP) + A(vP,vQ) + A(vQ,vP) + A(vQ,vQ); % Adding and grouping
            
            Ared = [A_II, A_IP; A_PI, A_PP];
        end
        
        function b_red = fullToReducedVector(obj,b)
            vF = obj.dof.free{1};
            vP = obj.dof.periodic_free;
            vQ = obj.dof.periodic_constrained;
            vI = setdiff(vF,vP);            
            b_I = b(vI);
            b_P = b(vP) + b(vQ);
            b_red = [b_I; b_P];
        end
        
        function b = reducedToFullVector(obj,bfree)
            vF = obj.dof.free;
            vP = obj.dof.periodic_free;
            vI = setdiff(vF{1},vP);
            
            b = zeros(obj.dof.ndof,1);
            b(vI) = bfree(1:1:size(vI,2));
            b(obj.dof.periodic_free) = bfree(size(vI,2)+1:1:size(bfree,1));
            b(obj.dof.periodic_constrained) = b(obj.dof.periodic_free);
            
        end
        
    end
    
end

