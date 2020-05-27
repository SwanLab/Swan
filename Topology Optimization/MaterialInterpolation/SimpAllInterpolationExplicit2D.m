classdef SimpAllInterpolationExplicit2D < SimpAllInterpolationExplicit
        
    methods  (Access = public)
        
        function obj = SimpAllInterpolationExplicit2D(cParams)
            obj.init(cParams);
            obj.nstre = 3;  
            obj.computeSymbolicInterpolationFunctions();
        end

    end    
  
end

