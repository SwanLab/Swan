classdef TopOpt_Problem_ComplianceLamPerimeter_st_Volume < TopOpt_Problem
    properties
    end
    methods
        function obj=TopOpt_Problem_ComplianceLamPerimeter_st_Volume(settings)
            obj@TopOpt_Problem(settings);
            obj.cost=ShFunc_CompliancePerimeter(settings);
            obj.constraint=ShFunc_Volume(settings.volume);             
        end
    end
end
