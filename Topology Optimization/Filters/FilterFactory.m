classdef FilterFactory < handle
    
    methods (Access = public)
        
        function obj = create(obj,type,optimizer)
            designVar = obj.getDesignVariableType(optimizer);
            switch type
                case 'P1'
                    switch designVar
                        case 'DENSITY'
                            obj = Filter_P1_Density();
                        case 'LEVELSET'
                            obj = Filter_P1_LevelSet();
                    end
                case 'PDE'
                    switch designVar
                        case 'DENSITY'
                            obj = Filter_PDE_Density();
                        case 'LEVELSET'
                            obj = Filter_PDE_LevelSet();
                    end
            end
        end
        
    end
    
    methods (Access = private, Static)
        
        function designVar = getDesignVariableType(optimizer)
            switch optimizer
                case {'MMA','PROJECTED GRADIENT','IPOPT'}
                    designVar = 'DENSITY';
                case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                    designVar = 'LEVELSET';
            end
        end
        
    end
    
end