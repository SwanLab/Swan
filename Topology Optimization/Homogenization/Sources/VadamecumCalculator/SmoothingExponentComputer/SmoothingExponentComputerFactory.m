classdef SmoothingExponentComputerFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type
                case 'Optimal'
                   obj = SmoothingExponentComputerOptimal(cParams); 
                case 'Given'
                   obj = SmoothingExponentComputerGiven(cParams); 
            end                    
        end
        
    end
    
end