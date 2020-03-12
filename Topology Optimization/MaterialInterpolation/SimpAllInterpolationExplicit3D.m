classdef SimpAllInterpolationExplicit3D < SimpAllInterpolationExplicit
        
    methods  (Access = public)
        
        function obj = SimpAllInterpolationExplicit3D(cParams)
            obj.init(cParams);
            obj.nstre = 6;  
            obj.computeSymbolicInterpolationFunctions();
        end

    end    
  
end

