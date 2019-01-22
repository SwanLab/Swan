classdef DesignVariableCreator < handle
    
    properties (Access = private)
        designVariable
        levelSet
    end
    
    methods (Access = public)
        
        function obj = DesignVariableCreator(settings,mesh)
            obj.createLevelSet(settings,mesh);
            obj.computeDesignVariable(settings.optimizer);
        end
        
        function x = getValue(obj)
            x = obj.designVariable;
        end
        
        function ls = getLevelSet(obj)
           ls = obj.levelSet;
        end
        
    end
    
    methods (Access = private)
       
        function createLevelSet(obj,settings,mesh)
            lsCreator    = LevelSetCreator.create(settings,mesh);
            obj.levelSet = lsCreator.getValue();
        end
        
        function computeDesignVariable(obj,optimizer)
            switch optimizer
                case {'SLERP','HAMILTON-JACOBI', 'PROJECTED SLERP'}
                    obj.designVariable = obj.levelSet;
                otherwise
                    phi = obj.levelSet;
                    obj.designVariable = 1 - heaviside(phi);
            end            
        end
        
    end
    
end

