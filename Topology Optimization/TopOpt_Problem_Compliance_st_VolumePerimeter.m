classdef TopOpt_Problem_Compliance_st_VolumePerimeter < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Compliance_st_VolumePerimeter(settings)
            obj@TopOpt_Problem(settings);
            obj.cost=ShFunc_Compliance;
            obj.constraint=ShFunc_VolumePerimeter(settings);
             
        end
    end
end
