classdef TopOpt_Problem < handle
    properties (GetAccess = public,SetAccess = public)
        cost
        constraint
        TOL
        x
        interpolation
        filter
        algorithm
        optimizer
        physicalProblem
        settings
    end
    methods (Access = public)
        function obj=TopOpt_Problem(settings)
            obj.settings=settings;
            obj.TOL=obj.settings.TOL;
            obj.physicalProblem=Physical_Problem(obj.settings.filename);        
            obj.interpolation=Interpolation.create(obj.TOL,settings.material,settings.method);                       
            switch obj.settings.optimizer
                case 'SLERP'
                    obj.optimizer=Optimizer_AugLag(settings,Optimizer_SLERP(settings));
                    obj.settings.ini_value=-0.7071;
                case 'PROJECTED GRADIENT'
                    obj.optimizer=Optimizer_AugLag(settings,Optimizer_PG(settings));
                    obj.settings.ini_value=1;
                case 'MMA'
                    obj.optimizer=Optimizer_MMA(settings);
                    obj.settings.ini_value=1;         
            end
            obj.filter=Filter.create(obj.settings.filter,obj.settings.optimizer);
        end
            
        function preProcess(obj)
            %initialize design variable
            obj.physicalProblem.preProcess;
            obj.filter.preProcess(obj.physicalProblem);
            switch obj.settings.initial_case
                case 'full'
                    obj.x=obj.settings.ini_value*ones(obj.physicalProblem.mesh.npnod,obj.physicalProblem.geometry.ngaus);
            end
        end
        function computeVariables(obj)
            obj.physicalProblem.computeVariables;
            obj.cost.computef(obj.x,obj.physicalProblem,obj.interpolation,obj.filter);
            obj.constraint.computef(obj.x, obj.physicalProblem, obj.interpolation,obj.filter);
            obj.x=obj.optimizer.solveProblem(obj.x,obj.cost,obj.constraint,obj.physicalProblem,obj.interpolation,obj.filter);
        end
        function postProcess(obj)
        end
    end
end
