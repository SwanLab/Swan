classdef Factory_BuilderDesignVarMonitor < handle
    
    methods (Access = public, Static)
        
        function builder = create(dim)
            switch dim
                case '2D'
                    builder = Builder_DesignVarMonitor_2D();
                case '3D'
                    builder = Builder_DesignVarMonitor_3D();
            end
        end
        
    end
    
end