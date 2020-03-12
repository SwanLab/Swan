classdef SimpallInterpolationExplicit2D < SimpallInterpolationExplicit
        
    methods  (Access = public)
        
        function obj = SimpallInterpolationExplicit2D(cParams)
            obj.init(cParams);
            obj.nstre = 3;  
            obj.ndim  = 2;   
            obj.computeSymbolicInterpolationFunctions();
        end

    end    
  
end

