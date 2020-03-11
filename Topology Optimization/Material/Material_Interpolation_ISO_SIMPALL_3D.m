classdef Material_Interpolation_ISO_SIMPALL_3D < Material_Interpolation_ISO_SIMPALL
    
    methods (Access = public)

        function obj= Material_Interpolation_ISO_SIMPALL_3D(cParams)
            obj.init(cParams)
        end
        
        function matProps=computeMatProp(obj, rho)
            dC=zeros(6,6,size(rho,1));
            if isempty(obj.mu_func)
                obj.computeSymProps();
                
                dmu = diff(obj.mu_sym);
                d = 3; % the dimension, 3 in 3D
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
            dC(1,3,:)= obj.dlam_func(rho);
            dC(2,1,:)= obj.dlam_func(rho);
            dC(2,2,:)= 2*obj.dmu_func(rho)+obj.dlam_func(rho);
            dC(2,3,:)= obj.dlam_func(rho);
            dC(3,1,:)= obj.dlam_func(rho);
            dC(3,2,:)= obj.dlam_func(rho);
            dC(3,3,:)= 2*obj.dmu_func(rho)+obj.dlam_func(rho);
            dC(4,4,:)= obj.dmu_func(rho);
            dC(5,5,:)= obj.dmu_func(rho);
            dC(6,6,:)= obj.dmu_func(rho);           
            
            matProps=struct;
            matProps.mu=mu;
            matProps.kappa=kappa;
            matProps.dC=dC;
        end
        
    end
end