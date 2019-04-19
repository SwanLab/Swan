classdef HomogenizedVarComputerFromVademecum ...
        < HomogenizedVarComputer        
    
    properties (Access = public)
        Ctensor
        density
    end
    
    
    methods (Access = public)
        
        function obj = HomogenizedVarComputerFromVademecum(cParams)
            s.fileName = cParams.fileName;
            v = VademecumVariablesLoader(s);
            v.load();
            obj.Ctensor = v.Ctensor;
            obj.density = v.density;            
        end        
        
        function computeMatProp(obj,rho)
            mx = sqrt(1-obj.rho);
            my = sqrt(1-obj.rho);
            Ctensor = obj.Ctensor.compute([mx,my]);
            
        end
        
    end
    
end