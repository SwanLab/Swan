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
                otherwise
                    error('Invalid interpolation order for element Tetrahedra.');
            end
        end
    end
end