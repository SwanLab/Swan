classdef Tetrahedra<Interpolation
    
    properties
    end
    
    methods
        % Constructor
        function obj = Tetrahedra(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'TETRAHEDRA';
            obj.order= 'LINEAR';
            obj.ndime = 3;
            obj.nnode = 4;
            obj.pos_nodes=[0 0 0;
                           1 0 0;
                           0 1 0;
                           0 0 1];
            obj.dvolu=1/6;
            obj.iteration=[1 1 1 2 2 3;
                   2 3 4 3 4 4];
        end
        function computeShapeDeriv(obj,posgp)
            obj.shape=[];
            s = posgp(1,:);
            t = posgp(2,:);
            u = posgp(3,:);
            obj.shape =[(ones(1,size(posgp,2))-t-s-u);
                s;
                t;
                u];
            obj.deriv=repmat([-1 1 0 0
                -1 0 1 0
                -1 0 0 1],1,1,size(posgp,2));            
        end
    end
    
end
