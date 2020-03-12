classdef SimpallInterpolationExplicit3D < SimpallInterpolationExplicit
        
    methods  (Access = public)
        
        function obj = SimpallInterpolationExplicit3D(cParams)
            obj.init(cParams);
            obj.nstre = 6;  
            obj.ndim  = 3;   
            obj.computeSymbolicInterpolationFunctions();
        end

    end    
  
end

