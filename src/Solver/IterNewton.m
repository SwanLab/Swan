classdef IterNewton < handle
    
    
    methods (Static)
        
        function xNew = solve(LHS,RHS,x)
            deltaX = -LHS\RHS;
            xNew = x + deltaX; 
        end
        
    end
    
end