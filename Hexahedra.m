classdef Hexahedra < Isoparametric
    properties
    end
    methods
        %constructor
        function obj=Hexahedra()
            obj = obj@Isoparametric();
            obj.type = 'HEXAHEDRA';
            obj.ndime = 3;          % 1D/2D/3D
            obj.nnode = 8;
            obj.ngaus = 8;          
            %COMPUTE WEIGP AND POSGP
      
            nlocs = 2;
            posgl(1)=-0.577350269189626;
            posgl(2)= 0.577350269189626;
            weigl(1)= 1.0;
            weigl(2)= 1.0;
            igaus=0;
            
            for ilocs=1:nlocs
                for jlocs=1:nlocs
                    for klocs=1:nlocs
                        igaus=igaus+1;
                        obj.weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs);
                        obj.posgp(1,igaus)=posgl(ilocs);
                        obj.posgp(2,igaus)=posgl(jlocs);
                        obj.posgp(3,igaus)=posgl(klocs);
                    end
                end
            end 
            for igauss=1:obj.ngaus
                s=obj.posgp(1,igauss);
                t=obj.posgp(2,igauss);
                u=obj.posgp(3,igauss);
                lcord(1,1)= -1; lcord(1,2)= -1; lcord(1,3)= -1;
                lcord(2,1)=  1; lcord(2,2)= -1; lcord(2,3)= -1;
                lcord(3,1)=  1; lcord(3,2)=  1; lcord(3,3)= -1;
                lcord(4,1)= -1; lcord(4,2)=  1; lcord(4,3)= -1;
                lcord(5,1)= -1; lcord(5,2)= -1; lcord(5,3)=  1;
                lcord(6,1)=  1; lcord(6,2)= -1; lcord(6,3)=  1;
                lcord(7,1)=  1; lcord(7,2)=  1; lcord(7,3)=  1;
                lcord(8,1)= -1; lcord(8,2)=  1; lcord(8,3)=  1;
                for inode=1:obj.nnode
                    obj.shape(inode,igauss)=(1+lcord(inode,1)*s)*(1+lcord(inode,2)*t)*(1+lcord(inode,3)*u)/8;
                    obj.deriv(1,inode,igauss)=lcord(inode,1)*(1+lcord(inode,2)*t)*(1+lcord(inode,3)*u)/8;
                    obj.deriv(2,inode,igauss)=lcord(inode,2)*(1+lcord(inode,1)*s)*(1+lcord(inode,3)*u)/8;
                    obj.deriv(3,inode,igauss)=lcord(inode,3)*(1+lcord(inode,1)*s)*(1+lcord(inode,2)*t)/8;
                end
            end
        end
    end

end