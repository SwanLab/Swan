classdef DesignVariable < handle & matlab.mixin.Copyable
    
    properties (GetAccess = public, SetAccess = protected)
        mesh
        type
    end
    
    properties (Access = public)
        valueOld                        
        value        
        nVariables
        rho
    end
    
    properties (Access = private)
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
                
        function restart(obj)
            obj.value = obj.valueOld;
        end
        
        function updateOld(obj)
            obj.valueOld = obj.value;
        end
        
        function objClone = clone(obj)
            objClone = copy(obj);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.type = cParams.type;
            obj.mesh = cParams.mesh;            
        end
        
        
    end
    
end

