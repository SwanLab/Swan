classdef SimpAllInterpolationImplicit2D < SimpAllInterpolationImplicit
        
    methods  (Access = public)
        
        function obj = SimpAllInterpolationImplicit2D(cParams)
            obj.init(cParams);
            obj.nstre = 3;  
            obj.computeSymbolicInterpolationFunctions();
        end

    end    
  
end