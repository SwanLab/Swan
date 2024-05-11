classdef Interpolator < handle

    methods (Static, Access = public)

        function interpolator = create(cParams)
            f = InterpolatorFactory();
            interpolator = f.create(cParams);
        end
    end    
end