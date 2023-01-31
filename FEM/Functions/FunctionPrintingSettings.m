classdef FunctionPrintingSettings < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        fun
    end
    
    methods (Access = public)

        function obj = FunctionPrintingSettings(cParams)
            obj.init(cParams);
        end
        
    end

    methods (Access = private)
        
        function init(obj,cParams)
            obj.fun = cParams.fun;
        end

    end
end

