classdef MaterialStokes < Material
    
    properties (Access = public)
        mu
        nElem
    end
    
    methods (Access = public) 
        
        function obj = MaterialStokes(cParams)
            obj.nElem = cParams.nelem;
        end

        function compute(obj)
            m = zeros(4,4,obj.nElem);
            m(1,1,:) = 1;
            m(2,2,:) = 1;
            m(3,3,:) = 1;
            m(4,4,:) = 1;
            obj.mu = m;
        end
        
    end
end