classdef NonLinearFilter < handle

    methods (Static, Access = public)
        function obj = create(cParams)
            f   = NonLinearFilterFactory();
            obj = f.create(cParams);
        end
    end

end