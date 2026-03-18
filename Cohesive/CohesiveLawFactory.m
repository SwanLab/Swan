classdef CohesiveLawFactory
    methods(Static)
        function law = create(cParams)
            switch cParams.lawType
                case 'TractionBiliniarUncoupled'
                    law = TractionBiliniarUncoupled(cParams);
            end
        end
    end
end
