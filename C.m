classdef C<Matrix_Local
    %Tensor containing the elastic constants of every element
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = C(nstre,nelem)
                material = Material_Elastic();
                epoiss = (material.matProp.kappa - material.matProp.mu)./(material.matProp.kappa + material.matProp.mu);
                eyoung = 4*material.matProp.kappa.*material.matProp.mu./(material.matProp.kappa + material.matProp.mu);
                obj.value = zeros(nstre,nstre,nelem);
                c1 = eyoung./(1-epoiss.^2);
                obj.value(1,1,:) = c1;
                obj.value(1,2,:) = c1.*epoiss;
                obj.value(2,1,:) = c1.*epoiss;
                obj.value(2,2,:) = c1;
                obj.value(3,3,:) = c1*0.5.*(1-epoiss);
        end
        
    end
    
end

