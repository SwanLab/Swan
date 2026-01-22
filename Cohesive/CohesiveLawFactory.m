classdef CohesiveLawFactory
    methods(Static)
        function law = create(lawType, params)
            switch lawType
                case 'Cubic'
                    law = CohesiveCubicLaw(params);
                case 'Bilinear'
                    law = CohesiveBilinearLaw(params);
                otherwise
                    error('Unknown law type %s', lawType);
            end
        end
    end
end
