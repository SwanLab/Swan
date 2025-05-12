classdef HardeningLawFactory < handle
    methods (Access = public, Static)
        function obj = create(cParams)
            switch cParams.type
                case 'Linear'
                    obj = HardeningLawLinear(cParams);
                case 'Exp'
                    obj = HardeningLawExp(cParams);
            end
        end 
    end
end