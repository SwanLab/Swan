classdef TopOpt_Problem_Compliance_st_Volume < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Compliance_st_Volume(settings)
            obj.TOL=settings.TOL;
            obj.cost=ShFunc_Compliance;
            obj.constraint=ShFunc_Volume(settings.Vfrac);
            obj.settings=settings;
            obj.filter=Filter_SLERP;
            %wip
            obj.physicalProblem=Physical_Problem(settings.filename);
            %wip
        end
    end
end
