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
            obj.cases(:,:,1)=[1 4 5;
                4 2 3;
                5 4 3];
            obj.cases(:,:,2)=[1 4 3;
                4 2 5;
                4 5 3];
            obj.cases(:,:,3)=[1 4 5;
                1 2 4;
                5 4 3];
            obj.selectcases =[1     0;
                2     0;
                3     3;
                0     2;
                0     1];
            obj.main_loop=[3 3];
            obj.extra_cases=[];
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
