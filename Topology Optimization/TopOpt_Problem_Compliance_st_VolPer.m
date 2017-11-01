classdef TopOpt_Problem_Compliance_st_VolPer < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Compliance_st_VolPer(x,TOL,physicalProblem,settings)
            obj.x=x;
            obj.TOL=TOL;
            obj.cost_func=ShFunc_Compliance;
            obj.constraint_func=ShFunc_VolumePerimeter(settings);
            obj.settings=settings;
            %wip
            obj.physicalProblem=physicalProblem;
            %wip
        end
    end
end
