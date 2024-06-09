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
                    q = Quadrature_Line();
                case 'TRIANGLE'
                    q = Quadrature_Triangle();
                case 'QUAD'
                    q = Quadrature_Quadrilateral();
                case 'TETRAHEDRA'
                    q = Quadrature_Tetrahedra();
                case 'HEXAHEDRA'
                    q = Quadrature_Hexahedra();
                otherwise
                    error('Invalid quadrature type.')
            end
            q.computeQuadrature(order);
        end
    end

    methods (Abstract, Access = public)
        computeQuadrature(obj,order)
    end

end
