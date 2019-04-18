classdef InactiveConstraintsModifierFactory < handle    
    
    methods (Access = public, Static)
        
        function obj = create(constraintCase)
            
            switch constraintCase
                case 'EQUALITY'
                    obj = InactiveConstraintsModifierNull();
                case 'INEQUALITY'
                    obj = InactiveConstraintsModifier();
                otherwise
                    error('Constraint case not valid.');
            end
            
        end
        
    end
               
end