classdef OrientationUpdater_MaxValue < OrientationUpdater
    
    methods (Access = protected)
        
        function computeOptimalIndexOrientation(obj)
            pS = obj.eigenValues;
            [~,ind] = max(pS);
            obj.optimalIndexOrientation = ind;
        end
        
    end
    
end