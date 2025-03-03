classdef OrientationUpdater_MaxAbsValue < OrientationUpdater
    
    methods (Access = protected)
        
        function computeOptimalIndexOrientation(obj)
            pS = obj.eigenValues;
            [~,ind] = max(abs(pS));
            obj.optimalIndexOrientation = ind;
        end
        
    end
    
end