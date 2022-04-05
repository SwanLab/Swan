classdef NewFilter < handle

    methods(Access = public, Static)

        function obj = create(cParams)
            f = NewFilterFactory();
            obj = f.create(cParams);
        end
        
    end

end