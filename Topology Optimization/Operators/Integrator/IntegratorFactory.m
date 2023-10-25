classdef IntegratorFactory < handle

    methods (Access = public, Static)

        function int = create(cParams)
            switch cParams.type
                case 'ScalarProduct'
                    int = IntegratorScalarProduct(cParams);
            end
        end

    end

end