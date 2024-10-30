classdef NonLinearFilterFactory < handle

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.type
                case 'Circle'
                    filter = NonLinearFilterCircle(cParams);
                case 'Ellipse'
                    %filter = NonLinearFilterEllipse(cParams);
                    filter = NonLinearFilterEllipsev2(cParams);
                case 'Segment'
                    filter = NonLinearFilterSegment(cParams);
                case 'Droplet'
                    filter = NonLinearFilterDroplet(cParams);
            end
        end

    end

end