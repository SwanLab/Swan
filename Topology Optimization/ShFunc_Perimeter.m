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
            obj.epsilon = 0.02;
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

        function computef(obj,x,physProblem,~,~)
            obj.checkFilterPre(physProblem);
            Msmooth=obj.filter_pde.Msmooth;
            obj.filter_pde.epsilon = obj.epsilon;
            x_reg=obj.filter_pde.getP1fromP1(x);
            rhs = obj.filter_pde.integrate_L2_function_with_shape_function(x);
            Perimeter = 0.5/obj.epsilon*((1 - x_reg)'*rhs);
            Perimeter_gradient = 0.5/obj.epsilon*(1 - 2*x_reg);
            
            constraint = Perimeter/obj.Perimeter_target - 1;
            constraint_gradient = Perimeter_gradient/obj.Perimeter_target;
            constraint_gradient = Msmooth*constraint_gradient;            
            
            obj.value=constraint;
            obj.gradient=constraint_gradient;
        end
        function checkFilterPre(obj, physicalProblem)
            if isempty(obj.filter_pde.Msmooth)
                nukn = 1;
                dof_per = DOF(physicalProblem.problemID,physicalProblem.geometry.nnode,physicalProblem.mesh.connec,nukn,physicalProblem.mesh.npnod,physicalProblem.mesh.scale);
                physicalProblem.dof = dof_per;
                switch physicalProblem.mesh.scale
                    case 'MACRO'
                        physicalProblem.dof.dirichlet = physicalProblem.dof.full_dirichlet;
                        physicalProblem.dof.dirichlet_values = physicalProblem.dof.full_dirichlet_values;
                        physicalProblem.dof.neumann = [];
                        physicalProblem.dof.neumann_values  = [];
                        physicalProblem.dof.constrained = physicalProblem.dof.compute_constrained_dof(physicalProblem.mesh.scale);
                        physicalProblem.dof.free = physicalProblem.dof.compute_free_dof();
                    case 'MICRO'
                        physicalProblem.dof.dirichlet = [];
                        physicalProblem.dof.dirichlet_values = [];
                        physicalProblem.dof.neumann = [];
                        physicalProblem.dof.neumann_values  = [];
                        physicalProblem.dof.constrained = physicalProblem.dof.compute_constrained_dof(physicalProblem.mesh.scale);
                        physicalProblem.dof.free = physicalProblem.dof.compute_free_dof();
                        
                        % Filter dofs coincides with physical_problem dofs
                end
                obj.filter_pde.preProcess(physicalProblem);
            end
        end
        
    end
end
