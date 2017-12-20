classdef ShFunc_Perimeter< Shape_Functional
    properties
        epsilon
        Perimeter_target
        filter_pde
    end
    methods
        function obj=ShFunc_Perimeter(settings)
            obj@Shape_Functional(settings);
            switch settings.perimeter.optimizer
                case 'SLERP'
                    obj.filter_pde=Filter_SLERP_PDE;
                otherwise
                    obj.filter_pde=Filter_Density_PDE;
            end
        end
        function Perimeter_target=get.Perimeter_target(obj)
            Perimeter_target=obj.target_parameters.Perimeter_target;
        end
        function computef(obj,x,physProblem,~,~)
            obj.checkFilterPre(physProblem)
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
        function epsilon=create_epsilon_perimeter(obj,coordinates,conectivities)
            xmin = min(coordinates);
            xmax = max(coordinates);
            epsilon0 = norm(xmax-xmin)/2;
            
            x1 = coordinates(conectivities(:,1));
            x2 = coordinates(conectivities(:,2));
            x3 = coordinates(conectivities(:,3));
            
            x1x2 = abs(x2-x1);
            x2x3 = abs(x3-x2);
            x1x3 = abs(x1-x3);
            hs = max([x1x2,x2x3,x1x3]');
            h = mean(hs);
            
            epsilon_end =h;
            frac = 2;
            kmax = ceil(log10(epsilon0/epsilon_end)/log10(frac));
            epsilon_iter = epsilon0./frac.^(1:kmax);
            epsilon=epsilon_iter(end);
        end
        function checkFilterPre(obj, physProblem)
            if isempty(obj.filter_pde.Msmooth)
                obj.filter_pde.preProcess(physProblem);
                obj.epsilon=obj.create_epsilon_perimeter(physProblem.mesh.coord,physProblem.mesh.connec);
            end
        end
        
    end
end
