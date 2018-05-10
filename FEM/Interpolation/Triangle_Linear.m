classdef Triangle_Linear<Interpolation 
    properties
    end
    methods
        % Constructor
        function obj = Triangle_Linear(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'TRIANGLE';
            obj.order = 'LINEAR';
            obj.ndime = 2;
            obj.nnode = 3;
            obj.pos_nodes = [0 0; 1 0; 0 1];
            obj.dvolu = 0.5;
        end
        function computeShapeDeriv(obj,posgp)
            obj.shape=[];
            for igaus=1:size(posgp,2)
                s = posgp(1,igaus);
                t = posgp(2,igaus);
                
                obj.shape(:,igaus) = [1.0-s-t;s;t];
                obj.deriv(:,:,igaus)=[-1.0 1.0 0.0;
                    -1.0 0.0 1.0];
            end
        end
    end
end
