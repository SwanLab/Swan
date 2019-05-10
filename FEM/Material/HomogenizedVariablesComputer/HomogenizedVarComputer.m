classdef HomogenizedVarComputer < handle
    
    properties (Access = public)
        C
        dC
        rho
        drho
    end
    
    properties (Access = protected)
       nelem
       designVariable
    end    
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = HomogenizedVarComputerFactory();
            obj = f.create(cParams);
        end
        
    end
        
end