classdef DensityCreator < handle
    
    properties (Access = protected)
        density
    end
    
    methods (Access = public)
                
        function d = getDensity(obj)
            d = obj.density;
        end
        
    end
    
end

