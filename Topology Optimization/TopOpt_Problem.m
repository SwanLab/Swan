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
        incremental_scheme
    end
    methods (Static)
        function obj=create(settings)
            switch settings.ptype
                case 'Compliance_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_Compliance_st_Volume(settings);
                case 'Compliance_st_VolumePerimeter'
                    settings.nconstr=2;
                    obj=TopOpt_Problem_Compliance_st_VolumePerimeter(settings);
                case 'ComplianceLamPerimeter_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_ComplianceLamPerimeter_st_Volume(settings);
                case 'Gripping'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_Gripping(settings);
                otherwise
                    disp('Problem not added')
            end
        end
    end
    methods (Access = public)
        function obj=TopOpt_Problem(settings)
            obj.settings=settings;
            obj.TOL=obj.settings.TOL;
            obj.physicalProblem=Physical_Problem(obj.settings.filename);    
            obj.incremental_scheme=Incremental_Scheme(obj.settings,obj.physicalProblem);
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
                case 'IPOPT'
                    obj.optimizer=Optimizer_IPOPT(settings);
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
                    obj.x=obj.settings.ini_value*ones(obj.physicalProblem.mesh.npnod,1);
            end
        end
        function computeVariables(obj)
            obj.physicalProblem.computeVariables;
            for istep = 1:obj.settings.nsteps
                disp(strcat('Incremental step: ', int2str(istep)))            
                obj.incremental_scheme.update_target_parameters(istep, obj.cost, obj.constraint, obj.optimizer);                              
                obj.cost.computef(obj.x,obj.physicalProblem,obj.interpolation,obj.filter);
                obj.constraint.computef(obj.x, obj.physicalProblem, obj.interpolation,obj.filter);
                obj.x=obj.optimizer.solveProblem(obj.x,obj.cost,obj.constraint,obj.physicalProblem,obj.interpolation,obj.filter);
            end
        end
        function postProcess(obj)
        end
    end
end
