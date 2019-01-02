classdef UnfittedMesh_Builder_Factory < handle
    methods (Access = public, Static)
        function concreteBuilder = create(type,ndim)
            switch type
                case 'INTERIOR'
                    switch ndim
                        case 2
                            concreteBuilder = UnfittedMesh_FlatSurface;
                        case 3
                            concreteBuilder = UnfittedMesh_Volumetric;
                    end
                case 'BOUNDARY'
                    switch ndim
                        case 2
                            concreteBuilder = UnfittedMesh_FlatCurve;
                        case 3
                            concreteBuilder = UnfittedMesh_3DSurface;
                    end
            end
        end
    end
end

