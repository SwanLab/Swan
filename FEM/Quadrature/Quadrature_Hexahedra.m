classdef Quadrature_Hexahedra<Quadrature
    properties
    end
    methods
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case 'LINEAR'
                    obj.ngaus = 8;
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
                otherwise
                    error('Invalid interpolation order for element Hexahedra.');
            end
        end
    end
end