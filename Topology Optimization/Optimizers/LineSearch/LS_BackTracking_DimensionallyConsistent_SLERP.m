classdef LS_BackTracking_DimensionallyConsistent_SLERP < LS_BackTracking_DimensionallyConsistent & LS_BackTracking_SLERP
    % NOTE: Besides this class is empty, fullfils a purpose: inheriting from 2 classes
    
    
    
    methods (Access = public)
        
        function obj = LS_BackTracking_DimensionallyConsistent_SLERP(cParams)
            obj.kfrac = cParams.kfrac;
        end
        
    end
end

