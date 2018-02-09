classdef TopOpt_Problem_Gripping < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Compliance_st_VolumePerimeter(settings)
            obj@TopOpt_Problem(settings);
            obj.cost=ShFunc_NotSelfAdjoint_Compliance(settings);
            obj.constraint=ShFunc_Volume(settings);
             
        end
    end
end
