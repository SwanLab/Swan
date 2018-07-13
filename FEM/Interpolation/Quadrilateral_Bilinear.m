classdef Quadrilateral_Bilinear < Interpolation
    properties
    end
    
    methods
        function obj = Quadrilateral_Bilinear(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'QUADRILATERAL';
            obj.order = 'LINEAR';
            obj.ndime = 2;
            obj.nnode = 4;
            obj.pos_nodes = [-1 -1; 1 -1; 1 1; -1 1];
            obj.dvolu = 4;
        end
        function computeShapeDeriv(obj,posgp)
            obj.shape=[];
            obj.deriv=[];
            s = posgp(1,:);
            t = posgp(2,:);
            obj.shape =[(ones(1,size(posgp,2))-t-s+s.*t)*0.25;
                (ones(1,size(posgp,2))-t+s-s.*t)*0.25;
                (ones(1,size(posgp,2))+t+s+s.*t)*0.25;
                (ones(1,size(posgp,2))+t-s-s.*t)*0.25];

            obj.deriv(1,1,:) =(-ones(1,size(posgp,2))+t)*0.25;
            obj.deriv(1,2,:) =(+ones(1,size(posgp,2))-t)*0.25;
            obj.deriv(1,3,:) =(+ones(1,size(posgp,2))+t)*0.25;
            obj.deriv(1,4,:) =(-ones(1,size(posgp,2))-t)*0.25;
            obj.deriv(2,1,:) =(-ones(1,size(posgp,2))+s)*0.25;
            obj.deriv(2,2,:) =(-ones(1,size(posgp,2))-s)*0.25;
            obj.deriv(2,3,:) =(+ones(1,size(posgp,2))+s)*0.25;
            obj.deriv(2,4,:) =(+ones(1,size(posgp,2))-s)*0.25;            
        end
    end
end

