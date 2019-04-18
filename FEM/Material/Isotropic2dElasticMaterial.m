classdef Isotropic2dElasticMaterial < IsotropicElasticMaterial

    properties (Access = protected)
    end
    
    methods (Access = public)
        
        function obj = Isotropic2dElasticMaterial(cParams)
            obj.nelem = cParams.nelem;
            obj.nstre = 3;            
            obj.createCtensor();            
        end
        
    end
    
    
    methods (Access = protected)
        
        function obj = computeC(obj)
            obj.computeYoungModulus();
            obj.computePoissonRatio();
            C = obj.C;            
            E = obj.E;
            nu = obj.nu;
            
            c1 = full(E./(1-nu.^2));
            C(1,1,:) = c1;
            C(1,2,:) = c1.*nu;
            C(2,1,:) = c1.*nu;
            C(2,2,:) = c1;
            C(3,3,:) = c1*0.5.*(1-nu);            
            obj.C = C;
        end
        
    end
    
    methods (Access = private)
        
        function computeYoungModulus(obj)
            k = obj.kappa;
            m = obj.mu;            
            E = 4*k.*m./(k + m);
            obj.E = E;
        end
        
        function computePoissonRatio(obj)
            k = obj.kappa;
            m = obj.mu;            
            nu = (k - m)./(k + m);
            obj.nu = nu;
        end
                
        
    end
    
end

