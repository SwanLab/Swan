classdef Solver < handle
    properties
        
    end
    
    methods (Access = ?Physical_Problem, Static)
        %Implement CREATE function when needed
        
        function stype = create()
            solver_type = 'DIRECT';
            switch solver_type
                case {'DIRECT'}
                    stype = Direct_solver();
                case {'ITERATIVE'}
                    error('Not implemented yet')
                otherwise
                    error('Invalid stype.')
            end
        end

    end
    
    


    
end
    
