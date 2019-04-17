classdef Material_Stokes < Material
    
    properties (Access = public)
        mu
    end
    
    methods (Access = public) 
        
        function obj = Material_Stokes(cParams)
            obj.nelem = cParams.nelem;
            mu = zeros(4,4,obj.nelem);
            mu(1,1,:) = 1;
            mu(2,2,:) = 1;
            mu(3,3,:) = 1;
            mu(4,4,:) = 1;
            obj.mu = mu;
        end
        
    end
end