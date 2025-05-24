classdef MaterialNavierStokes < Material
    
    properties (Access = public)
        nu
        nElem
    end
    
    methods (Access = public) 
        
        function obj = MaterialNavierStokes(cParams)
            obj.nElem   = cParams.nelem;
            obj.nuValue = cParams.nu;
            obj.compute();
        end

        function compute(obj)
            obj.nu = obj.nuValue * repmat(eye(4), [1, 1, obj.nElem]);
        end
 
    end
end