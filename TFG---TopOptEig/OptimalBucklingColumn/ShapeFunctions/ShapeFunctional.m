classdef ShapeFunctional < handle

    properties (Access = public)
        value
        gradient
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = ShapeFunctional_Factory();
            obj = f.create(cParams);
        end
        
    end
    
end