classdef IntegratorFactory < handle

    methods (Access = public, Static)

        function int = create(cParams)
            switch cParams.type
                case 'Function'
                    int = IntegratorFunction(cParams);
                case 'ScalarProduct'
                    int = IntegratorScalarProduct(cParams);
                case 'InternalEnergy'
                    int = IntegratorEnergy(cParams);
                case 'Error'
                    int = IntegratorError(cParams);
                case 'Unfitted'
                    int = IntegratorUnfitted(cParams);
            end
        end

    end

end