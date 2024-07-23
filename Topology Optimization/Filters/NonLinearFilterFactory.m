classdef NonLinearFilterFactory < handle

    methods (Access = public, Static)

        function filter = create(cParams)
            switch cParams.type
                case 'Circle'
                    filter = NonLinearFilterCircle(cParams);
                    % case Rectangle...
            end
        end

    end

end