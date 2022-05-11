classdef DimensionVariables < handle

    properties (Access = private)
    end

    methods (Static, Access = public)

        function obj = create(cParams)
            switch cParams.type
                case 'Scalar'
                    obj = DimensionScalar(cParams);
                case 'Vector'
                    obj = DimensionVector(cParams);
            end
        end
    end

end
