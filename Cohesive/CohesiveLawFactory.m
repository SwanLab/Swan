classdef CohesiveLawFactory
    methods(Static)
        function law = create(cParams)
            switch cParams.lawType
                case 'Cubic'
                    law = CohesiveCubicLaw(cParams);
                case 'Bilinear'
                    law = CohesiveBilinearLaw(cParams);
                case 'PPR'
                    law = CohesiveLawPPR(cParams);
                otherwise
                    error('Unknown law type %s', lawType);
            end
        end
    end
end
