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
            for igaus=1:size(posgp,2)
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                obj.shape(:,igaus) =[(1.-t-s+s*t)*0.25;
                    (1.-t+s-s*t)*0.25;
                    (1.+t+s+s*t)*0.25;
                    (1.+t-s-s*t)*0.25];
                
                obj.deriv(:,:,igaus) = [(-1.+t)*0.25 (+1.-t)*0.25 (+1.+t)*0.25 (-1.-t)*0.25;
                    (-1.+s)*0.25 (-1.-s)*0.25 (+1.+s)*0.25 (+1.-s)*0.25];
            end
        end
    end
end

