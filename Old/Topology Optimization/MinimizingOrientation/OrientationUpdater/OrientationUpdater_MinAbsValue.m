classdef OrientationUpdater_MinAbsValue < OrientationUpdater
    
    methods (Access = protected)
        
        function computeOptimalIndexOrientation(obj)
            pS = obj.eigenValues;
            [~,ind] = min(abs(pS));
            obj.optimalIndexOrientation = ind;
        end
        
    end
    
end