classdef Material_Interpolation_ISO_SIMPALL_2D < Material_Interpolation_ISO_SIMPALL
    properties        
    end
    methods
        function obj=Material_Interpolation_ISO_SIMPALL_2D(HSbounds)
            obj@Material_Interpolation_ISO_SIMPALL(HSbounds);
        end
        function matProps=computeMatProp(obj, rho)
            dC=zeros(3,3,size(rho,1));
            if isempty(obj.mu_func)
                obj.computeSymProps(rho);
                
                dmu = diff(obj.mu_sym);
                d = 2; % the dimension, 2 in 2D
                lam = obj.kappa_sym - 2/d*obj.mu_sym;
                dlam = diff(lam);
                
                obj.mu_func = matlabFunction(obj.mu_sym);
                obj.kappa_func = matlabFunction(obj.kappa_sym);
                obj.dmu_func=matlabFunction(dmu);
                obj.dlam_func=matlabFunction(dlam);
            end            
            mu= obj.mu_func(rho);
            kappa= obj.kappa_func(rho);
            dC(1,1,:)= 2*obj.dmu_func(rho)+obj.dlam_func(rho);
            dC(1,2,:)= obj.dlam_func(rho);
            dC(1,3,:)= zeros(length(rho),1);
            dC(2,1,:)= obj.dlam_func(rho);
            dC(2,2,:)= 2*obj.dmu_func(rho)+obj.dlam_func(rho);
            dC(2,3,:)= zeros(length(rho),1);
            dC(3,1,:)= zeros(length(rho),1);
            dC(3,2,:)= zeros(length(rho),1) ;
            dC(3,3,:)= obj.dmu_func(rho);            
            
            matProps=struct;
            matProps.mu=mu;
            matProps.kappa=kappa;
            matProps.dC=dC;
        end
        
    end
end