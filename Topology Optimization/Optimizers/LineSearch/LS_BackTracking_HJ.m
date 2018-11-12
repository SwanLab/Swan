classdef LS_BackTracking_HJ < LineSearch
    properties
        HJiter
        HJiter0
        HJiter_min = 1
    end
    
    methods
        function obj = LS_BackTracking_HJ(settings)
            obj.kappa = 1;
            obj.kappa_min = 1e-6;
            obj.kfrac = 2;
            obj.HJiter0 = settings.HJiter0;
            obj.HJiter = obj.HJiter0;
        end
    end
end

