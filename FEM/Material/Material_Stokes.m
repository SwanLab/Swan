classdef Material_Stokes < Material
    %Material_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mu
    end
    
    methods
        function obj = Material_Stokes(nelem)
            obj@Material(nelem);
            mu = zeros(4,4,obj.nelem);
            mu(1,1,:) = 1;
            mu(2,2,:) = 1;
            mu(3,3,:) = 1;
            mu(4,4,:) = 1;
            obj.mu = mu;
        end
    end
end