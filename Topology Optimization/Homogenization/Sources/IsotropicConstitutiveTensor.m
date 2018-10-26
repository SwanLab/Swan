classdef IsotropicConstitutiveTensor < FourthOrderTensor
    
    properties
        E
        nu
        kappa
        lambda
        mu
    end
    
    methods

        
    end
    
    methods (Access = public)
        
       function CreateIsotropicConstitutiveTensor(obj,E,nu)
            obj.E = E;
            obj.nu = nu;
            obj.computeLameParameters();
        end
        
        function computeLameParameters(obj)
            obj.computeShearModulus();
            obj.computeFirstLameParameter();
            obj.computeBulkModulus();
        end
        
        function computeShearModulus(obj)
            obj.mu = (obj.E/(2*(1+obj.nu)));
            
            if isUnit(obj.mu)
               obj.mu = simplify(obj.mu);
            end
                
        end
        
        function computeFirstLameParameter(obj)
            obj.lambda = (obj.E*obj.nu/((1+obj.nu)*(1-2*obj.nu)));
            
            if isUnit(obj.lambda)
               obj.lambda = simplify(obj.lambda);
            end
        end
        
        function Mu = getMu(obj)
            Mu = obj.mu;
        end
        
        function Lambda = getLambda(obj)
            Lambda = obj.lambda;
        end

    end
    
    methods (Abstract)
        computeBulkModulus(obj)
    end
    
    
    
end

