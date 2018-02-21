function [P_mas,P_menos] = polarization(igaus,element,problembsc,dim,chi,flag_change_micro)



if flag_change_micro == 1
    
    [Ce,Ce_plus,Ce_minus] = compute_consitutive_law(element,problembsc,igaus,dim,chi);
    gamma = element.material.opt_epsi;
    P_mas = zeros(dim.nstre,dim.nstre,dim.nelem);
    P_menos = zeros(dim.nstre,dim.nstre,dim.nelem);
    
    switch element.material.anisotropic
        
        case  'ANISOTROPIC'
            for ielem = 1:dim.nelem
                matprop.ce = Ce_plus(:,:,ielem);
                matprop.ci = Ce_minus(:,:,ielem);
                Polar_ani = polarization_ani(matprop);
                P_mas(:,:,ielem) = Polar_ani.Pe;
                P_menos(:,:,ielem) = Polar_ani.Pi;
            end
            
        case 'ISOTROPIC'
            mu = element.material.mu;
            lambda = element.material.lambda;
            [P_mas,P_menos] = polarization_iso(gamma,mu,lambda);
            
            for ielem = 1:dim.nelem
                P_mas(:,:,ielem) = P_mas;
                P_menos(:,:,ielem) = P_menos;
            end
    end
    
    
else
    P_mas = element.P_mas;
    P_menos = element.P_menos;
    
end

    
    
    


end