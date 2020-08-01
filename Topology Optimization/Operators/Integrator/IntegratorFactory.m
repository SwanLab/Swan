classdef IntegratorFactory
    methods (Access = public, Static)
        function integrator = create(cParams)
           switch cParams.type
                case 'COMPOSITE'
                    integrator = Integrator_Composite(cParams);
                case 'SIMPLE'
                    integrator = Integrator_Simple(cParams);
                case 'CutMesh'
                    integrator = IntegratorCutMesh(cParams);
               case 'Unfitted'
                    integrator = Integrator_Unfitted(cParams);
                otherwise
                    error('Invalid Unfitted Mesh type. Currently, integrator only works with INTERIOR and BOUNDARY.')
            end
        end
    end
end

