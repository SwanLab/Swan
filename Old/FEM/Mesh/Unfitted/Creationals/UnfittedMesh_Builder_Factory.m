classdef UnfittedMesh_Builder_Factory < handle
    
    methods (Access = public, Static)
        
        function concreteBuilder = create(type,ndim)
            switch ndim
                case 1
                    concreteBuilder = UnfittedMesh_StraightLine();
                case 2
                    switch type
                        case 'INTERIOR'
                            concreteBuilder = UnfittedMesh_FlatSurface();
                        case 'BOUNDARY'
                            concreteBuilder = UnfittedMesh_FlatCurve();
                    end
                case 3
                    switch type
                        case 'INTERIOR'
                            concreteBuilder = UnfittedMesh_Volumetric();
                        case 'BOUNDARY'
                            concreteBuilder = UnfittedMesh_3DSurface();
                    end
            end
        end
        
    end
    
end

