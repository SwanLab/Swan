classdef AbstractSettings < handle
    
    methods (Access = protected)
        
        function loadParams(obj,filename)
            cP = obj.getCustomParams(filename);
            obj.assignCustomParams(cP);
        end
        
    end
    
    methods (Access = private)
        
        function cP = getCustomParams(obj,filename)
            run(filename);
            vars = who;
            cP = struct;
            for i = 1:length(vars)
                propname = vars{i};
                if isprop(obj,propname)
                    cP.(propname) = eval(propname);
                else
%                     msg = [propname ' is not a property of ' class(obj) ' class.'];
%                     warning(msg);
                end
            end
        end
        
        function assignCustomParams(obj,cP)
            f = fields(cP);
            n = length(f);
            for i = 1:n
                param = f{i};
                obj.(param) = cP.(param);
            end
        end
        
    end
    
end
