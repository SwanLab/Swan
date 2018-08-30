classdef Material_Elastic_ISO_3D < Material_Elastic_ISO
    %Material_Elastic_ISO_3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Access = ?Material)
        function obj = Material_Elastic_ISO_3D(nelem)
            obj@Material_Elastic_ISO(nelem);
        end
    end
    
    methods (Access = protected)
        function obj = computeC(obj)
            C = zeros(6,6,obj.nelem);
            
            epoiss = (obj.kappa - obj.mu)./(obj.kappa + obj.mu);
            eyoung = 4*obj.kappa.*obj.mu./(obj.kappa + obj.mu);
            
            a = eyoung./((1+epoiss).*(1-2*epoiss));
            C(1,1,:) = a.*(1-epoiss);
            C(1,2,:) = a.*epoiss;
            C(2,1,:) = a.*epoiss;
            C(1,3,:) = a.*epoiss;
            C(3,1,:) = a.*epoiss;
            C(3,2,:) = a.*epoiss;
            C(2,3,:) = a.*epoiss;
            C(2,2,:) = a.*(1-epoiss);
            C(3,3,:) = a.*(1-epoiss);
            C(4,4,:) = a.*(1-2*epoiss)/2;
            C(5,5,:) = a.*(1-2*epoiss)/2;
            C(6,6,:) = a.*(1-2*epoiss)/2;
            
            obj.C = C;
        end
    end
    
end

