classdef Factory_BuilderDesignVarMonitor < handle
    
    methods (Access = public, Static)
        
        function builder = create(dim)
            switch dim
                case 1
                    builder = Builder_DesignVarMonitor_2D();
                case 2
                    builder = Builder_DesignVarMonitor_2D();
                case 3
                    builder = Builder_DesignVarMonitor_3D();
            end
        end
        
    end
    
end