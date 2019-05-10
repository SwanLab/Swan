classdef Isotropic3dElasticMaterial < IsotropicElasticMaterial
    
       methods (Access = public)
        
        function obj = Isotropic3dElasticMaterial(cParams)
            obj.init(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function obj = computeC(obj)
            obj.computeYoungModulus();
            obj.computePoissonRatio();
            C = obj.C;
            E = obj.E;
            nu = obj.nu;
            a = E./((1+nu).*(1-2*nu));
            Ch = a.*(1-nu);
            Cs = a.*nu;
            Css = a.*(1-2*nu)/2;
            C(1,1,:) = Ch;
            C(1,2,:) = Cs;
            C(2,1,:) = Cs;
            C(1,3,:) = Cs;
            C(3,1,:) = Cs;
            C(3,2,:) = Cs;
            C(2,3,:) = Cs;
            C(2,2,:) = Ch;
            C(3,3,:) = Ch;
            C(4,4,:) = Css;
            C(5,5,:) = Css;
            C(6,6,:) = Css;            
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

