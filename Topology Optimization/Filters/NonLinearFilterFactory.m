classdef NonLinearFilterFactory < handle

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.type
                case 'Circle'
                    filter = NonLinearFilterCircle(cParams);
                case 'Ellipse'
                    filter = NonLinearFilterEllipse(cParams);
            end
        end

    end

end