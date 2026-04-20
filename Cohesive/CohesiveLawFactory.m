classdef CohesiveLawFactory
    methods(Static)
        function law = create(cParams)
            switch cParams.lawType
                case 'TractionBiliniarUncoupled'
                    law = TractionBiliniarUncoupled(cParams);
                case 'TractionBiliniarCoupled'
                    law = TractionBiliniarCoupled(cParams);
            end
        end
    end
end
