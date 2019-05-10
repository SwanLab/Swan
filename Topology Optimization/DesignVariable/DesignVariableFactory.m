classdef DesignVariableFactory < handle
    
    properties (Access = private)
        cParams
    end
    
    methods (Access = public)
        
        function designVar = create(obj,cParams)
            obj.cParams = cParams;
            switch obj.cParams.type
                case 'LevelSet'                  
                    designVar = LevelSet(cParams);
                case 'Density'                  
                    designVar = Density(cParams);
                case 'MicroParams'
                    designVar = MicroParams(cParams);
            end
        end
        
    end    

end

