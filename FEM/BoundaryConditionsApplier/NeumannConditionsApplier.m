classdef NeumannConditionsApplier < BoundaryConditionsApplier
    
    
    methods (Access = public)

        
        function Ared = fullToReducedMatrix(obj,A)
            Ared = A;
        end
        
        function b_red = fullToReducedVector(obj,b)
            b_red = b;
        end
        
        function b = reducedToFullVector(obj,bfree)
           b = bfree;
        end
            
        
    end
  
    
end