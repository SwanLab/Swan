classdef OrientationUpdater_MinValue < OrientationUpdater
    
    methods (Access = protected)
        
        function computeOptimalIndexOrientation(obj)
            pS = obj.eigenValues;
            [~,ind] = min(pS);
            obj.optimalIndexOrientation = ind;
        end
        
    end
    
end