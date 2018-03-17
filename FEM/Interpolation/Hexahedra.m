classdef Hexahedra < Interpolation
    properties
    end
    methods
        %constructor
        function obj=Hexahedra(mesh)
            obj = obj@Interpolation(mesh);
            obj.type = 'HEXAHEDRA';
            obj.order = 'LINEAR';
            obj.ndime = 3;          % 1D/2D/3D
            obj.nnode = 8;
        end
        function computeShapeDeriv(obj,posgp)
            for igaus=1:size(posgp,2)                
                s=posgp(1,igaus);
                t=posgp(2,igaus);
                u=posgp(3,igaus);
                lcord(1,1)= -1; lcord(1,2)= -1; lcord(1,3)= -1;
                lcord(2,1)=  1; lcord(2,2)= -1; lcord(2,3)= -1;
                lcord(3,1)=  1; lcord(3,2)=  1; lcord(3,3)= -1;
                lcord(4,1)= -1; lcord(4,2)=  1; lcord(4,3)= -1;
                lcord(5,1)= -1; lcord(5,2)= -1; lcord(5,3)=  1;
                lcord(6,1)=  1; lcord(6,2)= -1; lcord(6,3)=  1;
                lcord(7,1)=  1; lcord(7,2)=  1; lcord(7,3)=  1;
                lcord(8,1)= -1; lcord(8,2)=  1; lcord(8,3)=  1;
                for inode=1:obj.nnode
                    obj.shape(inode,igaus)=(1+lcord(inode,1)*s)*(1+lcord(inode,2)*t)*(1+lcord(inode,3)*u)/8;
                    obj.deriv(1,inode,igaus)=lcord(inode,1)*(1+lcord(inode,2)*t)*(1+lcord(inode,3)*u)/8;
                    obj.deriv(2,inode,igaus)=lcord(inode,2)*(1+lcord(inode,1)*s)*(1+lcord(inode,3)*u)/8;
                    obj.deriv(3,inode,igaus)=lcord(inode,3)*(1+lcord(inode,1)*s)*(1+lcord(inode,2)*t)/8;
                end
            end
        end
    end
    
end