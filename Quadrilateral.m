classdef Quadrilateral < Isoparametric
    properties
    end
    methods
        %constructor
        function obj=Quadrilateral()
            obj = obj@Isoparametric();
            obj.type = 'QUAD';
            obj.ndime = 2;          % 1D/2D/3D
            obj.nnode = 4;
            obj.ngaus = 4;          % 
            %COMPUTE WEIGP AND POSGP
            nlocs = 2;
            posgl(1)=-0.577350269189626;
            posgl(2)= 0.577350269189626;
            weigl(1)= 1.0;
            weigl(2)= 1.0;
            igaus=0;
            for ilocs=1:nlocs
                for jlocs=1:nlocs
                    igaus=igaus+1;
                    obj.weigp(igaus)=weigl(ilocs)*weigl(jlocs);
                    obj.posgp(1,igaus)=posgl(ilocs);
                    obj.posgp(2,igaus)=posgl(jlocs);
                end
            end
            for igauss=1:obj.ngaus
            s=obj.posgp(1,igauss);
            t=obj.posgp(2,igauss);
            st=s*t;
            obj.shape(1,igauss)=(1.-t-s+st)*0.25;           %  4         3
            obj.shape(2,igauss)=(1.-t+s-st)*0.25;           %
            obj.shape(3,igauss)=(1.+t+s+st)*0.25;           %
            obj.shape(4,igauss)=(1.+t-s-st)*0.25;           %
            obj.deriv(1,1,igauss)=(-1.+t)*0.25;             %  1         2
            obj.deriv(1,2,igauss)=(+1.-t)*0.25;
            obj.deriv(1,3,igauss)=(+1.+t)*0.25;
            obj.deriv(1,4,igauss)=(-1.-t)*0.25;
            obj.deriv(2,1,igauss)=(-1.+s)*0.25;
            obj.deriv(2,2,igauss)=(-1.-s)*0.25;
            obj.deriv(2,3,igauss)=(+1.+s)*0.25;
            obj.deriv(2,4,igauss)=(+1.-s)*0.25; 
            end
        end
    end

end