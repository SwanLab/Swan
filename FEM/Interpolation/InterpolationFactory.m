classdef InterpolationFactory < handle
        
    methods (Access = public, Static)
        
        function obj = create(cParams)
            mesh  = cParams.mesh;
            order = cParams.order;
            switch mesh.type
                case 'EMPTY MESH'
                    obj = [];
                case 'LINE'
                    switch order
                        case 'LINEAR'
                            obj = Line_Linear(cParams);
                        otherwise
                            error('Invalid order for element LINE.');
                    end
                case 'TRIANGLE'
                    switch order
                        case 'LINEAR'
                            obj = Triangle_Linear(cParams);
                        case 'QUADRATIC'
                            obj = Triangle_Quadratic(cParams);
                        otherwise
                            error('Invalid order for element TRIANGLE.');
                    end
                case 'QUAD'
                    switch order
                        case 'LINEAR'
                            obj = Quadrilateral_Bilinear(cParams);
                        case 'QUADRATIC'
                            warning('PENDING TO BE TRASFORMED TO INTERPOLATION. SEE TRIANGLE_QUADRATIC AS EXAMPLE')
                            obj = Quadrilateral_Serendipity(cParamsr);
                        otherwise
                            error('Invalid order for element QUADRILATERAL.');
                    end
                case 'TETRAHEDRA'
                    obj = Tetrahedra_Linear(cParams);
                case 'HEXAHEDRA'
                    obj = Hexahedra_Linear(cParams);
                otherwise
                    error('Invalid mesh type.')
            end
        end
        
    end
    
end