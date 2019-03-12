classdef FilterFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.filter
                case 'P1'
                    switch cParams.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_P1_Density();
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            obj = Filter_P1_LevelSet();
                    end
                case 'PDE'
                    switch cParams.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_PDE_Density();
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            obj = Filter_PDE_LevelSet();
                    end
            end
        end
        
    end
    
end