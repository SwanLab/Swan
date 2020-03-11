classdef Material_Interpolation_ISO_SIMPALL_3D < Material_Interpolation_ISO_SIMPALL
    
    methods (Access = public)
        
        function obj= Material_Interpolation_ISO_SIMPALL_3D(cParams)
            obj.init(cParams)
            obj.nstre = 6;
            obj.ndim  = 3;
            obj.computeSymbolicInterpolationFunctions();
        end
        
    end
    
    
end
