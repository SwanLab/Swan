classdef ConstantLineSearch < handle
    
    properties (GetAccess = public)
        value
    end
    
    methods (Access = public)
        
        function obj = ConstantLineSearch(cParams)
            obj.value = cParams.value;
        end
        
    end
       
end