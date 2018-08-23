classdef Hexahedra_Linear < Interpolation
    properties
    end
    methods
        %constructor
        function obj=Hexahedra_Linear(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'HEXAHEDRA';
            obj.order = 'LINEAR';
            obj.ndime = 3;          % 1D/2D/3D
            obj.nnode = 8;
            obj.dvolu = 8;
            obj.pos_nodes=[-1 -1 -1;
                            +1 -1 -1;
                            +1 +1 -1;
                            -1 +1 -1;
                            -1 -1 +1;
                            +1 -1 +1;
                            +1 +1 +1;
                            -1 +1 +1];
            obj.iteration=[1 1 1 2 2 3 3 4 5 5 6 7;
                           2 4 5 3 6 4 7 8 6 8 7 8];
        end
        function computeShapeDeriv(obj,posgp)
            obj.shape=[];
            obj.deriv=[];            
            s=posgp(1,:);
            t=posgp(2,:);
            u=posgp(3,:);
            lcord(1,1)= -1; lcord(1,2)= -1; lcord(1,3)= -1;
            lcord(2,1)=  1; lcord(2,2)= -1; lcord(2,3)= -1;
            lcord(3,1)=  1; lcord(3,2)=  1; lcord(3,3)= -1;
            lcord(4,1)= -1; lcord(4,2)=  1; lcord(4,3)= -1;
            lcord(5,1)= -1; lcord(5,2)= -1; lcord(5,3)=  1;
            lcord(6,1)=  1; lcord(6,2)= -1; lcord(6,3)=  1;
            lcord(7,1)=  1; lcord(7,2)=  1; lcord(7,3)=  1;
            lcord(8,1)= -1; lcord(8,2)=  1; lcord(8,3)=  1;
            for inode=1:obj.nnode
                obj.shape(inode,:)=(ones(1,size(posgp,2))+lcord(inode,1)*s).*(ones(1,size(posgp,2))+lcord(inode,2)*t).*(ones(1,size(posgp,2))+lcord(inode,3)*u)/8;
                obj.deriv(1,inode,:)=lcord(inode,1).*(ones(1,size(posgp,2))+lcord(inode,2)*t).*(ones(1,size(posgp,2))+lcord(inode,3)*u)/8;
                obj.deriv(2,inode,:)=lcord(inode,2).*(ones(1,size(posgp,2))+lcord(inode,1)*s).*(ones(1,size(posgp,2))+lcord(inode,3)*u)/8;
                obj.deriv(3,inode,:)=lcord(inode,3).*(ones(1,size(posgp,2))+lcord(inode,1)*s).*(ones(1,size(posgp,2))+lcord(inode,2)*t)/8;
            end
        end
    end
    
end