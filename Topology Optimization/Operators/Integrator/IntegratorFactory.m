classdef IntegratorFactory
    methods (Access = public, Static)
        function integrator = create(mesh)
            switch mesh.unfittedType
                case 'INTERIOR'
                    integrator = Integrator_Interior;
                case 'BOUNDARY'
                    integrator = Integrator_Boundary;
                case 'COMPOSITE'
                    integrator = Integrator_Composite(mesh);
                otherwise
                    error('Invalid Unfitted Mesh type. Currently, integrator only works with INTERIOR and BOUNDARY.')
            end
        end
    end
end

