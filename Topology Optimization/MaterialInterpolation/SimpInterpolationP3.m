classdef SimpInterpolationP3 < SimpInterpolation
    
    methods (Access = public)
        
        function obj= SimpInterpolationP3(cParams)
            obj.init(cParams)
            obj.pExp = 3;
            obj.computeSymbolicInterpolationFunctions();
        end
           
    end
end