classdef MaterialNavierStokes < Material
    
    properties (Access = public)
        mu
        nElem
    end
    
    methods (Access = public) 
        
        function obj = MaterialNavierStokes(cParams)
            obj.nElem = cParams.nelem;
            obj.compute();
        end

        function compute(obj)
            nu = 1.51e-5;
            obj.mu = nu * repmat(eye(4), [1, 1, obj.nElem]);
        end
 
    end
end