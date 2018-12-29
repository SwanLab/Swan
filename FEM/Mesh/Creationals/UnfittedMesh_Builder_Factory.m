classdef UnfittedMesh_Builder_Factory < handle
    methods (Access = public, Static)
        function concreteBuilder = create(type,mesh_background,interpolation_background)
            ndim = mesh_background.ndim;
            switch type
                case 'INTERIOR'
                    switch ndim
                        case 2
                            concreteBuilder = UnfittedMesh_FlatSurface(mesh_background,interpolation_background);
                        case 3
                            concreteBuilder = UnfittedMesh_Volumetric(mesh_background,interpolation_background);
                    end
                case 'BOUNDARY'
                    switch ndim
                        case 2
                            concreteBuilder = UnfittedMesh_FlatCurve(mesh_background,interpolation_background);
                        case 3
                            concreteBuilder = UnfittedMesh_3DSurface(mesh_background,interpolation_background);
                    end
            end
        end
    end
end

