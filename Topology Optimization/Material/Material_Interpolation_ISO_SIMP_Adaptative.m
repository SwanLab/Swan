classdef Material_Interpolation_ISO_SIMP_Adaptative < Material_Interpolation_ISO_SIMP
    
    methods (Access = public)
        
        function obj= Material_Interpolation_ISO_SIMP_Adaptative(cParams)
            obj.init(cParams)
            obj.computeExponentP();
            obj.computeSymbolicInterpolationFunctions();
        end
       
    end
    
    methods (Access = private)
        
        function computeExponentP(obj)
           pUB = 2/(1-obj.nu1);
           pLB = 4/(1+obj.nu1);
           obj.pExp = max(pUB,pLB);
        end
    end
end