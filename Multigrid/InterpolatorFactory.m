classdef InterpolatorFactory < handle

    methods (Static)

        function obj = create(cParams)
            switch cParams.meshType
                case {'TRIANGLE'}
                    obj = TriangleInterpolator(cParams);
                case {'TETRAHEDRA'}
                    obj = TetrahedreInterpolator(cParams);
            end
        end
    end
end