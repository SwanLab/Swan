classdef Filter < handle

    methods(Access = public, Static)

        function obj = create(cParams)
            f = FilterFactory();
            obj = f.create(cParams);
        end

    end

end