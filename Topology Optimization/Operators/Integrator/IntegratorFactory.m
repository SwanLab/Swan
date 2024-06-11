classdef IntegratorFactory < handle

    methods (Access = public, Static)

        function int = create(cParams)
            switch cParams.type
                case 'Function'
                    int = IntegratorFunction(cParams);
                case 'Unfitted'
                    int = IntegratorUnfitted(cParams);
            end
        end

    end

end