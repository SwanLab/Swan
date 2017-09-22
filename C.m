classdef C<Matrix_Local
    %Tensor containing the elastic constants of every element
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = C(nstre,nelem, ndime)
                material = Material_Elastic();
                 obj.value = zeros(nstre,nstre,nelem);
                switch ndime
                    case 2
                        epoiss = (material.matProp.kappa - material.matProp.mu)./(material.matProp.kappa + material.matProp.mu);
                        eyoung = 4*material.matProp.kappa.*material.matProp.mu./(material.matProp.kappa + material.matProp.mu);
                       
                        c1 = eyoung./(1-epoiss.^2);
                        obj.value(1,1,:) = c1;
                        obj.value(1,2,:) = c1.*epoiss;
                        obj.value(2,1,:) = c1.*epoiss;
                        obj.value(2,2,:) = c1;
                        obj.value(3,3,:) = c1*0.5.*(1-epoiss);
                    case 3
                                a=material.matProp.kappa+2*material.matProp.mu;
                                obj.value(1,1,:) = a;
                                obj.value(1,2,:) = material.matProp.kappa;
                                obj.value(2,1,:) = material.matProp.kappa;
                                obj.value(1,3,:) = material.matProp.kappa;
                                obj.value(3,1,:) = material.matProp.kappa;
                                obj.value(3,2,:) = material.matProp.kappa;
                                obj.value(2,3,:) = material.matProp.kappa;
                                obj.value(2,2,:) = a;
                                obj.value(3,3,:) = a;
                                obj.value(4,4,:) = material.matProp.mu;
                                obj.value(5,5,:) = material.matProp.mu;
                                obj.value(6,6,:) = material.matProp.mu;
                end
        end
        
    end
    
end

