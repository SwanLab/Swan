classdef MomentumParameter < handle


    methods (Access = public, Static)

        function obj = create(cParams)
            f = MomentumParameterFactory();
            obj = f.create(cParams);
        end

    end

    methods (Access = public, Abstract)
        computeValue(obj)
    end

end