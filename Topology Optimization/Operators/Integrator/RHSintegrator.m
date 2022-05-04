classdef RHSintegrator < handle

    properties (Access = private)
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = RHSintegratorFactory();
            obj = f.create(s);
        end

    end
    
    methods (Access = public)
        
        function obj = RHSintegrator(cParams)
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
        end
        
    end
    
end