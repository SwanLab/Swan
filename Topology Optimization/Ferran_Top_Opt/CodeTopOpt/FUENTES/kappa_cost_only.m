function [kappa_opt,liter,vcostl,remesh] = kappa_cost_only(phifunct_l,g_nodal_n,theta_n,theta_t,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,alpha,h_C_0,Perim0,lambda_n,penalty_n,eta_n,iter,costfunc_n,vol_n,file_name,kappa_ini,Msmooth,Stiff_smooth,emass,flag_change_micro_ini,Group,file_gid,file_write,u0) 
            

freduc = element.material.kappa_reduction;
kappa_min = element.material.kappa_min;
kappa_end = element.material.kappa_end;

lmax = element.material.kappa_maxiter;
liter = 1; kappa_l = kappa_end; stop = 0;

vcostl = zeros(lmax,1);
remesh = 0;
flag_up_date_lambda = 0;
    while (liter <= lmax && stop==0)
             
            [phifunct_l1] = update_phifunc(theta_n,kappa_l,phifunct_l,g_nodal_n,norm_g_n,remesh,phifunct_l,problembsc.algorithm_update);
            [costk,thetak,~,volk] = equilibrium_update(phifunct_l1,alpha,h_C_0,Perim0,lambda_n,penalty_n,eta_n,theta_n,theta_t,element,problembsc,fixnodes,coordinates,fext,dim,Msmooth,Stiff_smooth,emass,flag_up_date_lambda,flag_change_micro_ini,Group,file_gid,file_write,u0);
                                                                                                                                                                          %    flag_up_date_lambda,flag_change_micro,Group,iter,file_gid,file_write,u_ant
            %[costk,thetak,volk] = call_new_cost(length(kappa)-ikappa,kappa(ikappa),theta_n,phifunct_l,g_nodal_n,norm_g_n,element,dim,problembsc,fixnodes,coordinates,fext,emass);
         %   cost_seg(ikappa,1) = costk;
          %  theta_seg(ikappa,1) = thetak;
          %  vol_seg(ikappa,1) = volk;
            cost_seg(liter) = costk;
            kappa(liter) = kappa_l;
            
        if (costk<costfunc_n)
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
        [~,i_min] = min(cost_seg);
        kappa_opt = kappa(i_min);
        %kappa_opt = kappa_min;
    end
end


    

