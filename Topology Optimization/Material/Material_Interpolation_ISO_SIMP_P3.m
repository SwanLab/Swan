classdef Material_Interpolation_ISO_SIMP_P3 < Material_Interpolation_ISO_SIMP
    
    methods (Access = public)
        
        function obj= Material_Interpolation_ISO_SIMP_P3(cParams)
            obj.init(cParams)
            obj.pExp = 3;
            obj.computeSymbolicInterpolationFunctions();
        end
           
    end
end