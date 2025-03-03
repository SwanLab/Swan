classdef LineSearchFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type
                case 'BACKTRACKING'
                    obj = LineSearch_Backtracking(cParams);          
                otherwise
                    error('%s is NOT a valid unconstrained optimizer.',cParams.optimizer);
            end
        end
        
    end
    
end