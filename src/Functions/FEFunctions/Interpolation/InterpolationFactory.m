classdef InterpolationFactory < handle

    methods (Access = public, Static)

        function obj = create(cParams)
            type  = cParams.type;
            order = cParams.order;
            switch type
                case 'EMPTY MESH'
                    obj = [];
                case 'LINE'
                    switch order
                        case 'CONSTANT'
                            obj = LineConstant(cParams);                        
                        case 'LINEAR'
                            obj = LineLinear(cParams);
                        case 'QUADRATIC'
                            obj = LineQuadratic(cParams);
                        case 'CUBIC'
                            obj = LineCubic(cParams);
                        otherwise
                            error('Invalid order for element LINE.');
                    end
                case 'TRIANGLE'
                    switch order
                        case 'CONSTANT' 
                            obj = TriangleConstant(cParams);
                        case 'LINEAR'
                            obj = TriangleLinear(cParams);
                        case 'QUADRATIC'
                            obj = TriangleQuadratic(cParams);
                        case 'CUBIC'
                            obj = TriangleCubic(cParams);
                        case 'RaviartThomas'
                            obj = Triangle_RaviartThomas(cParams);
                        case 'Nedelec'
                            obj = Triangle_Nedelec(cParams);
                        otherwise
                            error('Invalid order for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'CONSTANT'
                            obj = QuadrilateralConstant(cParams);
                        case 'LINEAR'
                            obj = QuadrilateralBilinear(cParams);
                        case 'QUADRATIC'
                            obj = QuadrilateralQuadratic(cParams);
                        case 'CUBIC'
                            obj = QuadrilateralCubic(cParams);
                        otherwise
                            error('Invalid order for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    switch order
                        case 'CONSTANT'
                            obj = TetrahedraConstant(cParams);
                        case 'LINEAR'
                            obj = TetrahedraLinear(cParams);
                        case 'QUADRATIC'
                            obj = TetrahedraQuadratic(cParams);
                        case 'CUBIC'
                            obj = TetrahedraCubic(cParams);
                        case 'RaviartThomas'
                            obj = Tetrahedra_RaviartThomas(cParams);
                        case 'Nedelec'
                            obj = Tetrahedra_Nedelec(cParams);
                        otherwise
                            error('Invalid order for element TETRAHEDRA.');
                    end
                case 'HEXAHEDRA'
                    switch order
                        case 'CONSTANT'
                            obj = HexahedraConstant(cParams);
                        case 'LINEAR'
                            obj = HexahedraLinear(cParams);
                        case 'QUADRATIC'
                            obj = HexahedraQuadratic(cParams);
                        case 'CUBIC'
                            obj = HexahedraCubic(cParams);
                        otherwise
                            error('Invalid order for element HEXAHEDRA.');
                    end
                otherwise
                    error('Invalid mesh type.')
            end
        end

    end

end