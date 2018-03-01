classdef ShFunc_Perimeter < Shape_Functional
    properties
        epsilon
        Perimeter_target
    end
    methods
        function obj = ShFunc_Perimeter(settings)
            if ~strcmp(settings.filter,'PDE')
                settings.filter = 'PDE';
                warning('Perimeter filter changed to PDE')
            end
            obj@Shape_Functional(settings);
            obj.Perimeter_target = settings.Perimeter_target;
            obj.target_parameters=settings.target_parameters;            
        end
        function epsilon=get.epsilon(obj)
            epsilon=obj.target_parameters.epsilon;
        end
        function computef(obj,x)
            obj.checkFilterPre(obj.filter);
            obj.filter.epsilon = obj.epsilon;
            x_reg = obj.filter.getP1fromP1(x);
            rhs = obj.filter.integrate_L2_function_with_shape_function(x);
            Perimeter = 0.5/obj.filter.epsilon*((1 - x_reg)'*rhs);
            Perimeter_gradient = 0.5/obj.filter.epsilon*(1 - 2*x_reg);
            
            constraint = Perimeter/obj.Perimeter_target - 1;
            constraint_gradient = Perimeter_gradient/obj.Perimeter_target;
            constraint_gradient = obj.filter.Msmooth*constraint_gradient;
            
            obj.value = constraint;
            obj.gradient = constraint_gradient;
        end
        %% !! DEBUGGING FUNCTION !! REMOVE WHEN DONE
        function checkFilterPre(obj, physicalProblem)
            if isempty(obj.filter.Msmooth)
                dof_phy = physicalProblem.dof;
                nukn = 1;
                dof_filter = DOF(physicalProblem.problemID,physicalProblem.geometry.nnode,physicalProblem.mesh.connec,nukn,physicalProblem.mesh.npnod,physicalProblem.mesh.scale);
                physicalProblem.dof = dof_filter;
                switch physicalProblem.mesh.scale
                    case 'MACRO'
                        dof_filter.dirichlet = physicalProblem.dof.full_dirichlet;
                        dof_filter.dirichlet_values = physicalProblem.dof.full_dirichlet_values;
                        dof_filter.neumann = [];
                        dof_filter.neumann_values  = [];
                        dof_filter.constrained = physicalProblem.dof.compute_constrained_dof(physicalProblem.mesh.scale);
                        dof_filter.free = physicalProblem.dof.compute_free_dof();
                    case 'MICRO'
                        dof_filter.dirichlet = [];
                        dof_filter.dirichlet_values = [];
                        dof_filter.neumann = [];
                        dof_filter.neumann_values  = [];
                        dof_filter.constrained = physicalProblem.dof.compute_constrained_dof(physicalProblem.mesh.scale);
                        dof_filter.free = physicalProblem.dof.compute_free_dof();
                        
                        % Filter dofs coincides with physical_problem dofs
                end
                physicalProblem.setDof(dof_filter)
                obj.filter.preProcess(physicalProblem);
                physicalProblem.setDof(dof_phy)
            end
        end
        
    end
end
