classdef Quadrature_Quadrilateral < Quadrature

    methods (Access = public)
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case 'CONSTANT'
                    obj.ngaus = 1;
                    obj.posgp(:,1) = [0,0];
                    obj.weigp = 4;
                    
                case 'LINEAR'
                    obj.ngaus = 4;
                    % Compute WEIGP and POSGP
                    a =  0.577350269189626;
                    obj.posgp(:,1) = [-a,-a];
                    obj.posgp(:,2) = [+a,-a];
                    obj.posgp(:,3) = [-a,+a];
                    obj.posgp(:,4) = [+a,+a];
                    
                    obj.weigp =  [1,1,1,1];%1*ones(1,obj.ngaus);

                case 'QUADRATICMSS' %SERENDIPITY, QUADRILATERAL QUADRATIC NOT IMPLEMENTED
                    obj.ngaus = 9;                    
                    % Compute WEIGP and POSGP
                    a =  0.77459667;
                    obj.posgp(:,1) = [ 0,+a];
                    obj.posgp(:,2) = [ 0, 0];
                    obj.posgp(:,3) = [+a,+a];
                    obj.posgp(:,4) = [-a,-a];
                    obj.posgp(:,5) = [-a, 0];
                    obj.posgp(:,6) = [+a, 0];
                    obj.posgp(:,7) = [+a,-a];
                    obj.posgp(:,8) = [-a,+a];
                    obj.posgp(:,9) = [ 0,-a];
                    
                    obj.weigp =(4/9)*ones(1,obj.ngaus);
                    
                case 'QUADRATICMASS'
                    posgl(1) =-0.774596669241483;
                    posgl(2) = 0.0;
                    posgl(3) = 0.774596669241483;
                    weigl(1) = 0.555555555555556;
                    weigl(2) = 0.888888888888889;
                    weigl(3) = 0.555555555555556;
                    obj.ngaus = 9;
                    igaus = 0;
                    nlocs = 3;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            igaus = igaus+1;
                            obj.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
                            obj.posgp(1,igaus) = posgl(ilocs);
                            obj.posgp(2,igaus) = posgl(jlocs);
                        end
                    end
                        
                otherwise
                    error('Invalid interpolation order for element Quadrilateral.');
            end
        end
    end
end