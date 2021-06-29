classdef SimpAllInterpolationImplicit3D < SimpAllInterpolationImplicit
    
    methods (Access = public)
        
        function obj= SimpAllInterpolationImplicit3D(cParams)
            obj.init(cParams)
            obj.nstre = 6;
            obj.computeSymbolicInterpolationFunctions();
        end
        
    end
    
    
end
