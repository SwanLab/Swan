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
            x_reg=obj.filter_pde.getP1fromP1(x,obj.epsilon);
            Perimeter = 0.5/obj.epsilon*((1 - x_reg)'*obj.filter_pde.rhs);
            Perimeter_gradient = 0.5/obj.epsilon*(1 - 2*x_reg);
            
            constraint = Perimeter/obj.Perimeter_target - 1;
            constraint_gradient = Perimeter_gradient/obj.Perimeter_target;
            constraint_gradient = Msmooth*constraint_gradient;            
            
            obj.value=constraint;
            obj.gradient=constraint_gradient;
        end
        function checkFilterPre(obj, physicalProblem)
            if isempty(obj.filter_pde.Msmooth)
                
                obj.dof_per=DOF(physicalProblem.problemID,physicalProblem.geometry.nnode,physicalProblem.mesh.connec,1,physicalProblem.mesh.npnod,physicalProblem.mesh.scale);
                
                switch physicalProblem.mesh.scale
                    case 'MACRO'
                        obj.dof_per.dirichlet = obj.dof_per.full_dirichlet;
                        obj.dof_per.dirichlet_values = obj.dof_per.full_dirichlet_values;
                    case 'MICRO'
                        % Filter dofs coincides with physical_problem dofs
                end

                obj.filter_pde.preProcess(physicalProblem,obj.dof_per);
            end
        end
        
    end
end
