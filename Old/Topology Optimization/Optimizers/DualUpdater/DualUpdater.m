classdef DualUpdater < handle
    
    properties (Access = public)
        dualVariable
    end
   
    methods (Access = public, Static)
    
        function obj = create(cParams)
            f   = DualUpdaterFactory();
            obj = f.create(cParams);
        end
    
    end
    
    methods (Access = public, Abstract)
        updateDualVariable(obj)
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.dualVariable = cParams.dualVariable;
        end
        
    end
    
    methods (Access = private)
        
        
        
    end
   
    
end