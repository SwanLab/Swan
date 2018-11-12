function [kappa_opt,vtheta] = kappa_bruteforce(element,dim,coordinates,problembsc,...
    phifunct_n,g_nodal_n,norm_g,fixnodes,fext,Msmooth,theta_n,iter,emass,vtheta)

% se busca el kappa que hace cost mas pequeño en un rango prefijado

liter = 1; 
kappa_ini = element.material.kappa_ini;
kappa_end = element.material.kappa_end;
lmax = element.material.kappa_maxiter;
hkappa = element.material.hkappa;
mkappa = size(hkappa,2);
if (iter<=mkappa+1)
    kappa_opt = hkappa(iter-1);
else
    while (liter <= lmax)
        kappa_l = (kappa_end-kappa_ini)*((liter-1)/(lmax-1))+kappa_ini;
        [phi_l] = update_phifunc(theta_n,kappa_l,phifunct_n,g_nodal_n,norm_g);
        
        % solve ku=f,topological derivative, and cost function
        [dummy1,dummy2,dummy3,dummy4,cost_l,theta_l,d_u,tstres,post,vol_omega,ener,vdisp,nbdata] = ...
            module_M(phi_l,element,fixnodes,problembsc,coordinates,fext,...
            dim,Msmooth,emass);
        
        vtheta(iter,liter,1) = liter;
        vtheta(iter,liter,2) = cost_l;
        vtheta(iter,liter,3) = theta_l*180/pi;
        vtheta(iter,liter,4) = kappa_l;
        vtheta(iter,liter,5) = vol_omega;
        vtheta(iter,liter,6) = ener.estra; % h_C
        vtheta(iter,liter,7) = ener.vol;   % L*|omega1|
        fprintf(1,'LITER: %3.0d COST_L %25.22f THETA_l %d KAPPA %25.19d %25.19d \n',liter,cost_l,theta_l*180/pi,kappa_l,vol_omega);
        liter=liter+1;
    end
    save('vtheta.mat','vtheta');
    kappa_opt=0.5;
end


    
end

