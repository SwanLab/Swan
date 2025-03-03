classdef HomogenizedVarComputer < handle
    
    properties (Access = public)
        Cref
        dCref
        C
        dC
        rho
        drho
        Pp
        dPp
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = HomogenizedVarComputerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public, Abstract)
        computeCtensor(obj)
        computePtensor(obj)
        computeDensity(obj)
    end
    
end