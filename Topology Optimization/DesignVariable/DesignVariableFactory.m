classdef DesignVariableFactory < handle
    
    properties (Access = private)
        cParams
        translator
        settingsDesignVar
    end
    
    methods (Access = public)
        
        function designVar = create(obj,cParams)            
            obj.cParams = cParams;
            obj.setupDesignVarSettings();
            switch obj.cParams.type
                case 'LevelSet'
                    designVar = LevelSet(obj.settingsDesignVar);
                case 'Density'
                    designVar = Density(obj.settingsDesignVar);
                case 'MicroParams'
                    designVar = MicroParams(obj.settingsDesignVar);
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function setupDesignVarSettings(obj)
            obj.settingsDesignVar.mesh  = obj.cParams.mesh;
            obj.settingsDesignVar.value = obj.cParams.value;
            obj.settingsDesignVar.type  = obj.cParams.type;
        end
        
    end
    
end

