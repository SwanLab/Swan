classdef LS_BackTracking_SLERP < LineSearch
    methods
        function obj = LS_BackTracking_SLERP
            obj.kappa = 1;
            obj.kappa_min = 1e-3;
            obj.kfrac = 2;
        end
    end
end

