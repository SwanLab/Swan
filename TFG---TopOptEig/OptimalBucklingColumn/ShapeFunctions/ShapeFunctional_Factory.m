classdef ShapeFunctional_Factory < handle

    properties (Access = private)
        designVariable
    end
    
    methods (Access = public)

        function sF = create(obj,cParams)
            obj.designVariable = cParams.designVariable;

            switch cParams.type
                case 'firstEignValue_functional'
                    sF = ShFunc_FirstEigenValue(cParams);
                case 'doubleEig1'
                    sF = Sh_doubleFirstEig(cParams);
                case 'doubleEig2'
                    sF = Sh_doubleSecondEig(cParams);
                case 'volume'
                    sF = Sh_volume(cParams);
                otherwise
                    error('Wrong cost name or not added to Cost Object');
            end
        end
    end
end
    