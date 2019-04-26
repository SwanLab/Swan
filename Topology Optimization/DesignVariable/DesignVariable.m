classdef DesignVariable < handle & matlab.mixin.Copyable
    
    properties (GetAccess = public, SetAccess = protected)
        meshGiD
        mesh        
        type
    end
    
    properties (Access = public)
        value        
    end
    
    methods (Access = public, Abstract)
        update(obj,value)
    end
    
    methods (Access = public, Static)
        
        function designVariable = create(cParams)
            f = DesignVariableFactory();
            designVariable = f.create(cParams);
        end        
        
    end
    
    methods (Access = public)
                
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type    = cParams.type;
            obj.meshGiD = cParams.mesh;            
        end
        
        
    end
    
end

