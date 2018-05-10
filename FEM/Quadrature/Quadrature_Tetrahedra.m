classdef Quadrature_Tetrahedra<Quadrature
    properties
    end
    methods
        function computeQuadrature(obj,order)
            switch order
                case 'LINEAR'
                    obj.ngaus = 1;          % tetrahedra
                    obj.weigp = 1/6;
                    obj.posgp = [1/4;1/4; 1/4];
                case 'QUADRATICMASS'
                    obj.ngaus=4;
                    a=0.58541020;
                    b=0.13819660;
                    obj.posgp=[a,b,b,b;
                                b,a,b,b;
                                a,a,b,b];
                    obj.weigp=[0.041666667,0.041666667,0.041666667,0.041666667];                          


                otherwise
                    error('Invalid interpolation order for element Tetrahedra.');
            end
        end
    end
end