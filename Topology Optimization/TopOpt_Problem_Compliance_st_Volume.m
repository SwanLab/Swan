classdef TopOpt_Problem_Compliance_st_Volume < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Compliance_st_Volume(settings)
            obj@TopOpt_Problem(settings);
            obj.cost=ShFunc_Compliance(settings);
            obj.constraint=ShFunc_Volume(settings);      
        end
    end
end
