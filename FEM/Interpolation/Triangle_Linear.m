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
            obj.deriv=[];            
            s = posgp(1,:);
            t = posgp(2,:);
            
            obj.shape = [ones(1,size(posgp,2))-s-t;s;t];
            obj.deriv=repmat([-1.0 1.0 0.0;
                -1.0 0.0 1.0],1,1,size(posgp,2));           
        end
    end
end
