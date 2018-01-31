classdef ShFunc_Perimeter< Shape_Functional
    properties
        epsilon
        Perimeter_target
        filter_pde
    end
    methods
        function obj=ShFunc_Perimeter(settings)
%            obj@Shape_Functional(settings);
            obj.Perimeter_target=settings.target_parameters.Perimeter_target;
            switch settings.perimeter.optimizer
                case 'SLERP'
                    obj.filter_pde=Filter_SLERP_PDE;
                otherwise
                    obj.filter_pde=Filter_Density_PDE;
            end
        end
%         function Perimeter_target=get.Perimeter_target(obj)
%             Perimeter_target=obj.target_parameters.Perimeter_target;
%         end
        function epsilon=get.epsilon(obj)
            epsilon=obj.target_parameters.epsilon;
        end
        function computef(obj,x,physProblem,~,~)
            obj.checkFilterPre(physProblem);
            Msmooth=obj.filter_pde.Msmooth;
            x_reg=obj.filter_pde.getP0fromP1_per(x,obj.epsilon);
            Perimeter = 0.5/obj.epsilon*((1 - x_reg)'*obj.filter_pde.rhs);
            Perimeter_gradient = 0.5/obj.epsilon*(1 - 2*x_reg);
            
            constraint = Perimeter/obj.Perimeter_target - 1;
            constraint_gradient = Perimeter_gradient/obj.Perimeter_target;
            constraint_gradient = Msmooth*constraint_gradient;            
            
            obj.value=constraint;
            obj.gradient=constraint_gradient;
        end
        function checkFilterPre(obj, physProblem)
            if isempty(obj.filter_pde.Msmooth)
                obj.filter_pde.preProcess(physProblem);
            end
        end
        
    end
end
