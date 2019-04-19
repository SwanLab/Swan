classdef HomegenizedVarComputerFromInterpolation ...
         < HomegenizedVarComputer
    
    properties (Access = public)
        matProps
    end
    
    properties (Access = private)
        interpolation
    end
    
    methods (Access = public)
        
        function obj = HomegenizedVarComputerFromInterpolation(cParams)
            int = Material_Interpolation.create(cParams);        
            obj.interpolation = int;
        end
        
        function mProps = computeMatProp(obj,rho)
           mProps = obj.interpolation.computeMatProp(rho);
           obj.matProps = mProps;
        end
        
        
    end
    
end