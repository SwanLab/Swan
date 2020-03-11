classdef Material_Interpolation < handle
    
    properties (SetAccess = protected, GetAccess = public)
        E_plus
        E_minus
        nu_plus
        nu_minus
        rho_plus
        rho_minus
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = MaterialInterpolationFactory;
            obj = f.create(cParams);
        end
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            cP = cParams.constitutiveProperties;            
            obj.rho_plus  = cP.rho_plus;
            obj.rho_minus = cP.rho_minus;
            obj.E_plus    = cP.E_plus;
            obj.E_minus   = cP.E_minus;
            obj.nu_plus   = cP.nu_plus;
            obj.nu_minus  = cP.nu_minus;            
        end
        
    end
    
end