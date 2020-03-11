classdef Material_Interpolation_ISO_SIMPALL_2D < Material_Interpolation_ISO_SIMPALL
        
    methods  (Access = public)
        
        function obj = Material_Interpolation_ISO_SIMPALL_2D(cParams)
            obj.init(cParams);
            obj.nstre = 3;  
            obj.ndim  = 2;   
            obj.computeSymbolicInterpolationFunctions();
        end

    end    
  
end