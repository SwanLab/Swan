function [kappa_opt,liter,vcostl] = kappa_cost_optim(phifunct_l,g_nodal_n,theta_n,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n)


kappa_min = element.material.kappa_min;
kappa_end = element.material.kappa_end;

lmax = element.material.kappa_maxiter;


fprintf(1,'ITER: %3.0d            COST_L %25.19d THETA  %10.6d  \n',iter,costfunc_n,theta_n*180/pi);


kappa = (linspace(kappa_min,kappa_end,100)).^1;
vcostl = size(kappa);
 for liter = 1:length(kappa) 
        kappa_l = kappa(liter);
               
        [phifunct_l1] = update_phifunc(theta_n,kappa_l,phifunct_l,g_nodal_n,norm_g_n);
        [phigp_l1] = interpol(phifunct_l1,element,dim,problembsc);
        
        
        
        
        
        [vol_l1] = cal_omega(phifunct_l1,dim,element,problembsc,coordinates);
        
        % solve ku=f
        [~,matCh] = module_M(phifunct_l1,phigp_l1,element,fixnodes,problembsc,coordinates,fext,dim);
        
        % cost function
        [cost_l1] = cal_cost_funct(matCh,problembsc,element,fext,vol_l1);
        vcostl(liter) = cost_l1;    
        
        figure(11)
        plot(kappa(1:liter),vcostl(1:liter),'-+')
        hold on
        plot([0 1],[costfunc_n costfunc_n],'r')
        hold off
        
        
        fprintf(1,'ITER: %3.0d LITER: %3.0d COST_L %25.19d THETA  %10.6d KAPPA %3.6f VOL %10.6d \n',...
                iter,liter,cost_l1,theta_n*180/pi,kappa_l,vol_l1);
 end
 
 [cost_l1,iter] = min(vcostl);
 kappa_opt = kappa(iter);
 
 if (cost_l1<costfunc_n)
   liter = 1;
 else
 liter=lmax;
 end
 
%kappa_opt = 0.3; 
 
end


    

