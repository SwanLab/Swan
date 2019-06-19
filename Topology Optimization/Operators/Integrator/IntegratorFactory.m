classdef IntegratorFactory
    methods (Access = public, Static)
        function integrator = create(cParams)
            switch cParams.mesh.unfittedType
                case 'INTERIOR'
                    integrator = Integrator_Interior(cParams);
                case 'BOUNDARY'
                    integrator = Integrator_Boundary(cParams);
                case 'COMPOSITE'
                    integrator = Integrator_Composite(cParams);
                case 'SIMPLE'
                    integrator = Integrator_Simple(cParams);
                otherwise
                    error('Invalid Unfitted Mesh type. Currently, integrator only works with INTERIOR and BOUNDARY.')
            end
        end
    end
end

