classdef Material_Elastic_ISO_3D < Material_Elastic_ISO
    %Material_Elastic_ISO_3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function obj = Material_Elastic_ISO_3D(nelem)
            obj = obj@Material_Elastic_ISO(nelem);
            obj = obj.computeC;
        end
        
%         function obj = setProps(obj,props)
%             obj = setProps_parent(obj,props);
%             obj = obj.computeC;
%         end
%     end
%     
%     methods (Access = protected)
%         function obj = setProps_parent(obj,props)
%             setProps_parent@Material_Elastic_ISO(obj,props);
%         end
     end
    
    methods (Access = private)
        function obj = computeC(obj)
            C = zeros(6,6,obj.nelem);
            
            a = obj.kappa+2*obj.mu;
            C(1,1,:) = a;
            C(1,2,:) = obj.kappa;
            C(2,1,:) = obj.kappa;
            C(1,3,:) = obj.kappa;
            C(3,1,:) = obj.kappa;
            C(3,2,:) = obj.kappa;
            C(2,3,:) = obj.kappa;
            C(2,2,:) = a;
            C(3,3,:) = a;
            C(4,4,:) = obj.mu;
            C(5,5,:) = obj.mu;
            C(6,6,:) = obj.mu;
            
            obj.C = C;
        end
    end
    
end

