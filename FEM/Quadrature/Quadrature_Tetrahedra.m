classdef Quadrature_Tetrahedra < Quadrature
    methods
        function computeQuadrature(obj,order)
            computeQuadrature@Quadrature(obj,order);
            switch order
                case 'CONSTANT'
                    obj.ngaus = 1;          % tetrahedra
                    obj.weigp = 1/6;
                    obj.posgp = [1/4;1/4;1/4];                    
                    
                case 'LINEAR'
                    obj.ngaus = 1;          % tetrahedra
                    obj.weigp = 1/6;
                    obj.posgp = [1/4;1/4;1/4];
                    
                case 'QUADRATIC'
                    obj.ngaus = 4;
                    a = 0.58541020;
                    b = 0.13819660;
                    obj.posgp = [a,b,b,b;
                        b,a,b,b;
                        b,b,a,b];
                    obj.weigp = [0.041666667,0.041666667,0.041666667,0.041666667];
                    
                case 'CUBIC'
                    obj.ngaus = 5;
                    a = 0.25;
                    b = 0.5;
                    c = 1/6;
                    obj.posgp = [a,b,c,c,c;
                        a,c,c,c,b;
                        a,c,c,b,c];
                    obj.weigp = 1/6*[-0.8,0.45,0.45,0.45,0.45];
                    
                case 'QUADRATICMASS'
                    obj.ngaus = 4;
                    a = 0.58541020;
                    b = 0.13819660;
                    obj.posgp = [a,b,b,b;
                        b,a,b,b;
                        a,a,b,b];
                    obj.weigp = [0.041666667,0.041666667,0.041666667,0.041666667];
                    
                otherwise
                    error('Invalid interpolation order for element Tetrahedra.');
            end
        end
    end
end