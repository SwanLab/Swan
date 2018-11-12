function [kappa_opt,liter,vcostl] = kappa_cost_only_old(phifunct_l,g_nodal_n,theta_n,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n)

freduc = element.material.kappa_reduction;

kappa_min = element.material.kappa_min;
kappa_end = element.material.kappa_end;

lmax = element.material.kappa_maxiter;
liter = 1; kappa_l = kappa_end; stop = 0;

fprintf(1,'ITER: %3.0d            COST_L %25.19d THETA  %10.6d  \n',iter,costfunc_n,theta_n*180/pi);
    while (liter <= lmax && stop==0)
    
               
        [phifunct_l1] = update_phifunc(theta_n,kappa_l,phifunct_l,g_nodal_n,norm_g_n);
        [phigp_l1] = interpol(phifunct_l1,element,dim,problembsc);
        
        [vol_l1] = cal_omega(phifunct_l1,dim,element,problembsc,coordinates);
        
        % solve ku=f
        [~,matCh] = module_M(phifunct_l1,phigp_l1,element,fixnodes,problembsc,coordinates,fext,dim);
        
        % cost function
        [cost_l1] = cal_cost_funct(matCh,problembsc,element,fext,vol_l1);
        vcostl(liter) = cost_l1;     
        
        fprintf(1,'ITER: %3.0d LITER: %3.0d COST_L %25.19d THETA  %10.6d KAPPA %3.6f VOL %10.6d \n',...
                iter,liter,cost_l1,theta_n*180/pi,kappa_l,vol_l1);
            
        if (cost_l1<costfunc_n)
            kappa_opt = kappa_l;
            stop = 1;
        else
            kappa_l = freduc*kappa_l;
            if (kappa_l < kappa_min)
                kappa_l = kappa_min;
                liter=lmax;
            end
%             fprintf(1,'ITER: %3.0d LITER: %3.0d COST_L %25.22f THETA_l %d KAPPA %25.19d VOL %25.19d \n',...
%                 iter,liter,cost_l,theta_l*180/pi,kappa_l,vol_l);
            liter=liter+1;
        end      
        %constr = 1 - vol_l1/element.material.Vfrac;
        %element = update_lambda(problembsc,element,constr);
    end
    if (liter > lmax) % maximo numero de iteraciones
        kappa_opt = kappa_min;
    end
end


    

