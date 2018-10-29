classdef IsotropicConstitutiveTensor2D < IsotropicConstitutiveTensor
    
    properties
    end
    
    methods
        
        function obj = IsotropicConstitutiveTensor2D(E,nu)
            obj = obj@IsotropicConstitutiveTensor(E,nu);
        end
        
        function computeBulkModulus(obj)
            obj.kappa = obj.lambda + obj.mu;
        end
        

        
    end
    
end

