classdef SimpInterpolationAdaptative < SimpInterpolation
    
    methods (Access = public)
        
        function obj= SimpInterpolationAdaptative(cParams)
            obj.init(cParams)
            obj.computeExponentP();
            obj.computeSymbolicInterpolationFunctions();
        end
       
    end
    
    methods (Access = private)
        
        function computeExponentP(obj)
           pUB = 2/(1-obj.matProp.nu1);
           pLB = 4/(1+obj.matProp.nu1);
           obj.pExp = max(pUB,pLB);
        end
    end
end