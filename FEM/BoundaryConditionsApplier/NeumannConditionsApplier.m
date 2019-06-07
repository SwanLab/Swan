classdef NeumannConditionsApplier < handle
    
    
    methods (Access = public)

        
        function Ared = full_matrix_2_reduced_matrix(obj,A)
            Ared = A;
        end
        
        function b_red = full_vector_2_reduced_vector(obj,b)
            b_red = b;
        end
        
        function b = reduced_vector_2_full_vector(obj,bfree)
           b = bfree;
        end
            
        
    end
  
    
end