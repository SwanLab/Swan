classdef Optimizer < handle
   
    properties
        has_converged = false;
        stop_vars
        target_parameters = struct;
        nconstr
        constraint_case
    end
    
    
    methods
        function obj = Optimizer(settings)
            obj.nconstr   = settings.nconstr;
            obj.target_parameters = settings.target_parameters;
            obj.constraint_case=settings.constraint_case;
        end
    end
    
    
    methods (Access = protected)
        function cons = setConstraint_case(obj,constraint)
            cons = constraint;
            switch obj.constraint_case
                case 'EQUALITY'
                case 'INEQUALITY'
                    contr_inactive_value = -obj.objfunc.lambda(:)./obj.objfunc.penalty(:);
                    inactive_constr = contr_inactive_value' > constraint.value;
                    cons.value(inactive_constr) = contr_inactive_value(inactive_constr);
                    cons.gradient(:,inactive_constr) = 0;
                otherwise
                    error('Constraint case not valid.');
            end
        end
        
        
    end
end