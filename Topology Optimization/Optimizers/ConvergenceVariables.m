classdef ConvergenceVariables < handle
    
    properties (GetAccess = public, SetAccess = private)
        values
        nVar
    end
    
    methods (Access = public)
        
        function obj = ConvergenceVariables(n)
            obj.nVar = n;
        end
        
        function reset(obj)
            obj.values = [];
            obj.nVar = 0;
        end
        
        function append(obj,value)
            obj.values(obj.nVar+1) = value;
            obj.nVar = obj.nVar + 1;
        end
        
        function set(obj,value)
            obj.values(end) = value;
        end
        
    end
    
end