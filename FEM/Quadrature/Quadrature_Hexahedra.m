classdef Quadrature_Hexahedra < Quadrature
    
    methods
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case {'LINEAR','CONSTANT'}
                    obj.ngaus = 8;
                    nlocs = 2;
                    posgl(1) = -0.577350269189626;
                    posgl(2) = 0.577350269189626;
                    weigl(1) = 1.0;
                    weigl(2) = 1.0;
                    igaus = 0;
                    
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            for klocs = 1:nlocs
                                igaus = igaus+1;
                                obj.weigp( igaus) = weigl(ilocs)*weigl(jlocs)*weigl(klocs);
                                obj.posgp(1,igaus) = posgl(ilocs);
                                obj.posgp(2,igaus) = posgl(jlocs);
                                obj.posgp(3,igaus) = posgl(klocs);
                            end
                        end
                    end
                    
                case 'QUADRATICMASS'
                    posgl(1) = -0.774596669241483;
                    posgl(2) = 0.0;
                    posgl(3) = 0.774596669241483;
                    weigl(1) = 0.555555555555556;
                    weigl(2) = 0.888888888888889;
                    weigl(3) = 0.555555555555556;
                    obj.ngaus = 27;
                    igaus = 0;
                    nlocs = 3;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            for klocs = 1:nlocs
                                igaus = igaus+1;
                                obj.weigp( igaus) = weigl(ilocs)*weigl(jlocs)*weigl(klocs);
                                obj.posgp(1,igaus) = posgl(ilocs);
                                obj.posgp(2,igaus) = posgl(jlocs);
                                obj.posgp(3,igaus) = posgl(klocs);
                            end
                        end
                    end
                    
                otherwise
                    error('Invalid interpolation order for element Hexahedra.');
            end
        end
    end
end