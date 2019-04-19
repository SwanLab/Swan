classdef Material < handle
    
    properties (Access = protected)
        nelem
    end
    
    
    methods (Access = public, Static)
        function material = create(cParams)
            f = MaterialFactory();
            material = f.create(cParams);
        end
    end
    
end