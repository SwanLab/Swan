classdef TopOpt_Problem < handle
    properties
        cost
        constraint
        cost_gradient
        constraint_gradient
        cost_func
        constraint_func
        TOL
        x
        method
        physicalProblem
        settings
    end
    methods
        function preProcess(obj)
            %choose interpolation
            switch obj.settings.material
                case 'ISOTROPIC'
                    switch obj.settings.method
                        case 'SIMPALL'
                            obj.method=Interpolation_ISO_SIMPALL(obj.TOL);
                        otherwise
                            disp('Method not added')
                    end
            end        
        end
        function computeVariables(obj)
            obj.preProcess;
            [obj.cost, obj.cost_gradient]=obj.cost_func.computef(obj.x,obj.physicalProblem,obj.method);
            [obj.constraint, obj.constraint_gradient]=obj.constraint_func.computef(obj.x,obj.physicalProblem,obj.interpolation);
        end
        postProcess(obj)
    end
end