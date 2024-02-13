classdef Quadrature_Hexahedra < Quadrature
    
    methods
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case 'CONSTANT'
                    obj.ngaus = 1;
                    obj.weigp = 1;
                    obj.posgp = [0;0;0];
                case 'LINEAR'
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
                    
                case 'QUADRATIC'
                    % Copied from QUADRATICMASS
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
                    
                case 'ORDER10'
                    posgl = [0.973906528517172;0.865063366688985;0.679409568299024;0.433395394129247;0.148874338981631;-0.148874338981631;-0.433395394129247;-0.679409568299024;-0.865063366688985;-0.973906528517172];
                    weigl = [0.0666713443086878;0.149451349150581;0.219086362515982;0.269266719309996;0.295524224714753;0.295524224714753;0.269266719309996;0.219086362515982;0.149451349150581;0.0666713443086878];
                    obj.ngaus = 1000;
                    igaus = 0;
                    nlocs = 10;
                    for ilocs = 1:nlocs
                        for jlocs = 1:nlocs
                            for klocs = 1:nlocs
                                igaus = igaus+1;
                                obj.weigp(igaus) = weigl(ilocs)*weigl(jlocs)*weigl(klocs);
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