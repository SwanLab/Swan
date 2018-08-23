classdef Solver < handle
    properties
        
    end
    
    methods (Static)
        function stype = create()
            solver_type = 'DIRECT';
            switch solver_type
                case {'DIRECT'}
                    stype = Direct_Solver();
                case {'CHOLEVSKY'}
                    % At least up to ndof ~5e4, Direct is still faster
                    stype = Cholesky_Direct_Solver();
                case {'ITERATIVE'}
                    error('Not implemented yet')
                otherwise
                    error('Invalid stype.')
            end
        end
    end
end
    
