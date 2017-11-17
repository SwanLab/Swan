classdef TopOpt_Problem_Compliance_st_Volume < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Compliance_st_Volume(settings)
            obj.TOL=settings.TOL;
            obj.cost_func=ShFunc_Compliance;
            obj.constraint_func=ShFunc_Volume(settings);
            obj.settings=settings;
            %wip
            obj.physicalProblem=Physical_Problem(settings.filename);
            %wip
        end
    end
end
