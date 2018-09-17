classdef Material_Elastic_ISO_2D < Material_Elastic_ISO
    %Material_Elastic_ISO_3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods 
        function obj = Material_Elastic_ISO_2D(nelem)
            obj@Material_Elastic_ISO(nelem);
        end
    end
    
    methods (Access = protected)
        function obj = computeC(obj)
            C = zeros(3,3,obj.nelem);
            
            epoiss = (obj.kappa - obj.mu)./(obj.kappa + obj.mu);
            eyoung = 4*obj.kappa.*obj.mu./(obj.kappa + obj.mu);
            
            c1 = full(eyoung./(1-epoiss.^2));
            C(1,1,:) = c1;
            C(1,2,:) = c1.*epoiss;
            C(2,1,:) = c1.*epoiss;
            C(2,2,:) = c1;
            C(3,3,:) = c1*0.5.*(1-epoiss);
            
            obj.C = C;
        end
    end
    
end

