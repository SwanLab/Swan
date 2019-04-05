classdef DesignVariableFactory < handle
    
    properties (Access = private)
        cParams
        translator = OptimizerToDesignVariableTranslator();
        settingsDesignVar
    end
    
    methods (Access = public)
        
        function designVar = create(obj,cParams)
            obj.cParams = cParams;
            obj.setupDesignVarSettings();
            switch obj.getDesignVarType()
                case 'LevelSet'
                    designVar = LevelSet(obj.settingsDesignVar);
                case 'Density'
                    designVar = Density(obj.settingsDesignVar);
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function setupDesignVarSettings(obj)
            obj.settingsDesignVar.mesh = obj.cParams.mesh;
            obj.settingsDesignVar.value = obj.cParams.value;
        end
        
        function type = getDesignVarType(obj)
            type = obj.translator.translate(obj.cParams.optimizer);
        end
        
    end
    
end

