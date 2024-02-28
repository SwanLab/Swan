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
                        case 'LINEAR'
                            obj = Line_Linear(cParams);
                        case 'CONSTANT'
                            obj = Line_Constant(cParams);
                        otherwise
                            error('Invalid order for element LINE.');
                    end
                case 'TRIANGLE'
                    switch order
                        case 'CONSTANT' 
                            obj = Triangle_Constant(cParams);
                        case 'LINEAR'
                            obj = Triangle_Linear(cParams);
                        case 'QUADRATIC'
                            obj = Triangle_Quadratic(cParams);
                        case 'CUBIC'
                            obj = Triangle_Cubic(cParams);
                        otherwise
                            error('Invalid order for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'CONSTANT'
                            obj = Quadrilateral_Constant(cParams);
                        case 'LINEAR'
                            obj = Quadrilateral_Bilinear(cParams);
                        case 'QUADRATIC'
                            obj = Quadrilateral_Quadratic(cParams);
%                             obj = Quadrilateral_Serendipity(cParams);
                        case 'CUBIC'
                            obj = Quadrilateral_Cubic(cParams);
                        otherwise
                            error('Invalid order for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    switch order
                        case 'CONSTANT'
                            obj = Tetrahedra_Constant(cParams);
                        case 'LINEAR'
                            obj = Tetrahedra_Linear(cParams);
                        case 'QUADRATIC'
                            obj = Tetrahedra_Quadratic(cParams);
                        case 'CUBIC'
                            obj = Tetrahedra_Cubic(cParams);
                        otherwise
                            error('Invalid order for element TETRAHEDRA.');
                    end
                case 'HEXAHEDRA'
                    switch order
                        case 'CONSTANT'
                            obj = Hexahedra_Constant(cParams);
                        case 'LINEAR'
                            obj = Hexahedra_Linear(cParams);
                        case 'QUADRATIC'
                            obj = Hexahedra_Quadratic(cParams);
                        case 'CUBIC'
                            obj = Hexahedra_Cubic(cParams);
                        otherwise
                            error('Invalid order for element HEXAHEDRA.');
                    end
                otherwise
                    error('Invalid mesh type.')
            end
        end

    end

end