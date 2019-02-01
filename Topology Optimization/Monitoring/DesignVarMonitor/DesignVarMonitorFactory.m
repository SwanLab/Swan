classdef DesignVarMonitorFactory < handle
    
    methods (Access = public)
        
        function monitor = create(obj,shallDisplay,settings,mesh)
            if shallDisplay
                switch obj.designVariable(settings.optimizer)
                    case 'Density'
                        switch settings.pdim
                            case '2D'
                                monitor = DesignVarMonitor_Density_2D(mesh);
                            case '3D'
                                monitor = DesignVarMonitor_Density_3D(mesh);
                        end
                        
                    case 'LevelSet'
                        
                    otherwise
                        error('Invalid Design Variable')
                end
            else
                monitor = DesignVarMonitor_Null(mesh);
            end
            
        end
        
    end
    
    methods (Access = private, Static)
        
        function var = designVariable(optimizer)
            switch optimizer
                case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                    var = 'LevelSet';
                case {'PROJECTED GRADIENT','MMA','IPOPT'}
                    var = 'Density';
            end
            
        end
        
    end
    
end