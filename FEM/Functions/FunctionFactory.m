classdef FunctionFactory < handle

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.functionType
                case 'P0'
                    obj = P0Function(cParams);
                case 'P1'
                    obj = P1Function(cParams);
                case 'P1Disc'
                    obj = P1DiscontinuousFunction(cParams);
                case 'FGaus'
                    obj = FGaussDiscontinuousFunction(cParams);
            end
        end

    end
end