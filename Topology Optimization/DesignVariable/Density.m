classdef Density < DesignVariable
    
    properties (Access = private)
        levelSetCreatorSettings
    end
    
    methods (Access = public)
        
        function obj = Density(cParams)
            obj.init(cParams);
            obj.levelSetCreatorSettings = cParams.levelSetCreatorSettings;
            obj.createValue();
        end
        
        function update(obj,value)
            obj.value = value;
        end
        
    end
    
    methods (Access = private)
        
        function createValue(obj)
            s = obj.levelSetCreatorSettings;
            lsCreator  = LevelSetCreator.create(s);
            phi        = lsCreator.getValue();
            obj.value  = 1 - heaviside(phi);
        end
        
    end
    
end

