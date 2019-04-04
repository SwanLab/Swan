classdef AbstractSettings < handle
    
    properties (Access = protected, Abstract)
        defaultParamsName
    end
    
    properties (Access = private)
        customParams
    end
    
    methods (Access = protected)
        
        function obj = AbstractSettings()
            obj.loadParams(obj.defaultParamsName);
        end
        
        function loadParams(obj,param_filenme)
            run(param_filenme);
            obj.customParams = who;
            obj.clearCustomParams();
            
            for i = 1:length(obj.customParams)
                param = obj.customParams{i};
                if isprop(obj,param)
                    obj.(param) = eval(param);
                else
                    obj.warnOfInvalidCustomParams(param);
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function clearCustomParams(obj)
            obj.removeVar('param_filenme');
            obj.removeVar('obj');
        end
        
        function warnOfInvalidCustomParams(obj,param)
            warning([param ' is not a property of ' class(obj)]);
        end
        
        function removeVar(obj,var)
            pos = strcmp(obj.customParams,var);
            obj.customParams(pos) = [];
        end
        
    end
    
end