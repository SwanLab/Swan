classdef HardeningLawFactory < handle
    methods (Access = public, Static)
        function obj = create(cParams)
            switch cParams.type
                case 'Linear'
                    obj = HardeningLawLinear(cParams);
                case 'Exponential'
                    obj = HardeningLawExp(cParams);
                case 'AT2'
                    obj = HardeningLawAT2(cParams);
                case 'AT1'
                    obj = HardeningLawAT1(cParams);
            end
        end
    end
end