classdef DesignVariableTruss < handle
    
    properties (Access = public)
        value
        valueOld
    end
    
    properties (Access = private)
        ubRad
        lbRad
        ubThick
        lbThick
    end
    
    methods (Access = public)
        
        function obj = DesignVariableTruss(cParams)
            obj.ubRad   = cParams.ubR;
            obj.lbRad   = cParams.lbR;
            obj.ubThick = cParams.ubT;
            obj.lbThick = cParams.ubT;
        end
        
        function init(obj,x0)
            obj.value = x0;
        end
        
    end
    
    methods (Access = public)
        
        function update(obj,x)
            varN = length(x)/2;
            xR   = x(1:varN);
            xT   = x(varN+1:end);
            xR   = min(obj.ubRad,max(obj.lbRad,xR));
            xT   = min(obj.ubThick,max(obj.lbThick,xT));
            obj.value = [xR;xT];
        end
        
        function updateOld(obj)
            obj.valueOld = obj.value;
        end
       
    end
    
end