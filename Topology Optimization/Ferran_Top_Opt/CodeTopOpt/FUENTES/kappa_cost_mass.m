function [kappa_opt] = kappa_cost_mass(element,dim,coordinates,problembsc,...
    phifunct_n,g_nodal_n,norm_g,fixnodes,fext,Msmooth,theta_n,iter,emass,costfunc_n,vol_n)

% se busca el kappa que cumple los siguientes criterios:
% Dado el incremento se fija el intervalo de busqueda, por ejemplo:
% - iter=1  --> 0.5
% - iter=2  --> 0.5
% - iter>=3 --> 0.05


kappa_min = element.material.kappa_min;
ndata = size(element.material.kappa_end,2);
if (iter-1<=ndata)
    kappa_end = element.material.kappa_end(iter-1);
else
    kappa_end = element.material.kappa_end(ndata);
end
freduc = element.material.kappa_reduction;
ndata = size(element.material.vol_reduction,2);
if (iter-1<=ndata)
    vol_reduc = element.material.vol_reduction(iter-1);
else
    vol_reduc = element.material.vol_reduction(ndata);
end
lmax = element.material.kappa_maxiter;
liter = 1; kappa_l = kappa_end; stop = 0;
    while (liter <= lmax && stop==0)
        [phi_l] = update_phifunc(theta_n,kappa_l,phifunct_n,g_nodal_n,norm_g);
        % solve ku=f,topological derivative, and cost function
        [dummy1,dummy2,dummy3,dummy4,cost_l,theta_l,d_u,tstres,post,vol_l,ener,vdisp,nbdata] = ...
            module_M(phi_l,element,fixnodes,problembsc,coordinates,fext,...
            dim,Msmooth,emass);
        
        vtheta(liter,1) = liter;
        vtheta(liter,2) = cost_l;
        vtheta(liter,3) = theta_l*180/pi;
        vtheta(liter,4) = kappa_l;
        vtheta(liter,5) = vol_l;
        vtheta(liter,6) = ener.estra; % h_C
        vtheta(liter,7) = ener.vol;   % L*|omega1|
        %fprintf(1,'LITER: %3.0d COST_L %25.22f THETA_l %d KAPPA %25.19d %25.19d \n',liter,cost_l,theta_l*180/pi,kappa_l,vol_l);
    
        cond_cost=0; cond_vol=0; cond_theta=0;
        if (cost_l<costfunc_n)
            cond_cost = 1;
        end
        vol_min = vol_n - vol_reduc*vol_n/100;
        if (vol_min< vol_l && vol_l<=vol_n)
            cond_vol = 1;
        end
        if (theta_l<theta_n)
            cond_theta = 1;
        end
        
        if (cond_cost==1 && cond_vol==1 && cond_theta==1)
            kappa_opt = kappa_l;
            stop = 1;
        elseif (cond_cost==1 && cond_vol==1)
            kappa_opt = kappa_l;
            stop = 1;
        elseif (cond_cost==1 && cond_theta==1)
            kappa_opt = kappa_l;
            stop = 1;
        elseif (cond_vol==1)
            kappa_opt = kappa_l;
            stop = 1;
        else
            kappa_l = freduc*kappa_l;
            if (kappa_l<kappa_min)
                kappa_l = kappa_min;
                liter=lmax;
            end
            liter=liter+1;
        end
            
    end
    if (liter > lmax) % maximo numero de iteraciones
        kappa_opt = kappa_min;
    end
        
    save('vtheta.mat','vtheta');
    
end


    

