classdef IntegratorFactory
    methods (Access = public, Static)
        function integrator = create(mesh)
            switch mesh.meshType
                case 'INTERIOR'
                    integrator = Integrator_Interior;
                case 'BOUNDARY'
                    integrator = Integrator_Boundary;
                case 'COMPOSITE'
                    integrator = Integrator_Composite(mesh);
                otherwise
                    error('Invalid Mesh type. Currently, integrator only works with INTERIOR and BOUNDARY.')
            end
        end
    end
end

