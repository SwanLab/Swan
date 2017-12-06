classdef ShFunc_Perimeter< Shape_Functional
    properties 
        epsilon
        Perimeter_target
    end
    methods
        function obj=ShFunc_Perimeter(persettings)
            obj.epsilon=persettings.epsilon;
            obj.Perimeter_target=persettings.Perimeter_target;
        end
    end
    methods
        function computef(obj,x,~,~,filter)
            
            %             switch algorithm
            %                 case 'level_set'
            %                     rhs = faireF2(coordinates',element.conectivities',design_variable);
            %
            %                 otherwise % case 'Projected_gradient'
            %                     rhs = Msmooth*design_variable;
            %
            %             end
            Msmooth=filter.Msmooth;
            x_reg=filter.getP0fromP1_perimeter(x,obj.epsilon);
            Perimeter = 0.5/obj.epsilon*((1 - x_reg)'*filter.rhs);
            Perimeter_gradient = 0.5/obj.epsilon*(1 - 2*x_reg);
            
            constraint = Perimeter/obj.Perimeter_target - 1;
            constraint_gradient = Perimeter_gradient/obj.Perimeter_target;
            constraint_gradient = Msmooth*constraint_gradient;
            
            
            obj.value=constraint;
            obj.gradient=constraint_gradient;
        end
    end
end
