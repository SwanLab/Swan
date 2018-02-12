classdef Solver < handle
    properties
        
    end
    
    methods (Access = ?Physical_Problem, Static)
        %Implement CREATE function when needed
        
        function stype = create(ptype)
            switch ptype
                case {'ELASTIC','THERMAL'}
                    stype = Solver_Dirichlet_Conditions();
                case 'HYPERELASTIC'
                    stype = Solver_NewtonRaphson();
                otherwise
                    error('Invalid stype.')
            end
        end

    end
    
end
    
