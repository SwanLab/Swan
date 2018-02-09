classdef TopOpt_Problem_Gripping < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_Gripping(settings)
            obj@TopOpt_Problem(settings);
            obj.cost=ShFunc_NonSelfAdjoint_Compliance(settings);
            obj.constraint=ShFunc_Volume(settings);
             
        end
    end
end
