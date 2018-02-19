classdef Quad_Mass<Isoparametric
    %Triangle_Linear Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Constructor
        function obj = Quad_Mass
            obj = obj@Isoparametric;
            obj.type = 'QUAD';
            obj.ndime = 2;
            obj.nnode=4;
            obj.ngaus=9;
            
            posgl(1)=-0.774596669241483;
            posgl(2)= 0.0;
            posgl(3)= 0.774596669241483;
            weigl(1)= 0.555555555555556;
            weigl(2)= 0.888888888888889;
            weigl(3)= 0.555555555555556;
            
            igaus=0;
            nlocs=3;
            for ilocs=1:nlocs
                for jlocs=1:nlocs
                    igaus=igaus+1;
                    obj.weigp(  igaus)=weigl(ilocs)*weigl(jlocs);
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