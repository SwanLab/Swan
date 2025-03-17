classdef InactiveConstraintsModifierFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            
            switch cParams.type
                case 'EQUALITY'
                    obj = InactiveConstraintsModifierNull(cParams);
                case 'INEQUALITY'
                    obj = InactiveConstraintsModifier(cParams);
                otherwise
                    error('Constraint case not valid.');
            end
            
        end
        
    end
    
end