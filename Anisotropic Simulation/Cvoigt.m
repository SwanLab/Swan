classdef Cvoigt < handle
    
    methods (Access = public, Static)

        function tensor = create(cParams)
            f = TensorFactory();
            tensor = f.create(cParams);
        end

    end

end