classdef Quadrature < handle

    properties (GetAccess = public, SetAccess = protected)
        posgp
        weigp
        ngaus
        order
    end
    
    methods (Static, Access = public)

        function q = create(mesh, order)
            switch mesh.type
                case 'LINE'
                    q = QuadratureLine();
                case 'TRIANGLE'
                    q = QuadratureTriangle();
                case 'QUAD'
                    q = QuadratureQuadrilateral();
                case 'TETRAHEDRA'
                    q = QuadratureTetrahedra();
                case 'HEXAHEDRA'
                    q = QuadratureHexahedra();
                otherwise
                    error('Invalid quadrature type.')
            end
            q.computeQuadrature(order);
        end
    end

    methods (Abstract, Access = protected)
        computeQuadrature(obj,order)
    end

end
