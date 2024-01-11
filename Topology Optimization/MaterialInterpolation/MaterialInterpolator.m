classdef MaterialInterpolator < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f   = MaterialInterpolatorFactory;
            obj = f.create(cParams);
        end
        
    end    
    
end