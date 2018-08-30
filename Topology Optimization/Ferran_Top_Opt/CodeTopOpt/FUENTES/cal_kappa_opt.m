function [ kappa_opt,liter,vcostl,remesh,iter_sub] = cal_kappa_opt(phifunct_n,g_nodal_n,theta_n,theta_t,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,weights,h_C_0,Perim0,lambda_n,penalty_n,eta_n,iter,costfunc_n,vol_n,file_name,kappa_ini,Msmooth,Stiff_smooth,emass,iter_sub,d_u,flag_change_micro_ini,Group,file_gid,file_write,u0)
                                                     
% Different functions to obtain kappa_opt
switch problembsc.meth_kappa_opt
    case 'IMPOSED'
        [kappa_opt] = kappa_imposed(element,iter);
    case 'BRUTE_FORCE'
        [kappa_opt,vtheta] = kappa_bruteforce(element,dim,coordinates,problembsc,...
            phifunct_n,g_nodal_n,norm_g,fixnodes,fext,Msmooth,theta_n,iter,emass,vtheta);
    case 'FREQUENCY'
        [kappa_opt] = kappa_frequency(dim,phifunct_n,g_nodal_n,norm_g,theta_n,...
            element,coordinates,coordinates,problembsc);            
    case 'COST_MASS_THETA'
        [kappa_opt] = kappa_cost_mass(element,dim,coordinates,problembsc,...
            phifunct_n,g_nodal_n,norm_g,fixnodes,fext,Msmooth,theta_n,iter,emass,costfunc_n,vol_n);
    case 'COST_ONLY'
        [kappa_opt,liter,vcostl,remesh] = kappa_cost_only(phifunct_n,g_nodal_n,theta_n,theta_t,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n);
  
    case 'COST_ONLY_REMESH'
        [kappa_opt,liter,vcostl,remesh,iter_sub] = kappa_cost_only_remesh(phifunct_n,g_nodal_n,theta_n,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n,vol_n,file_name,kappa_ini,emass,iter_sub);

    case 'COST_OPT'
        %[kappa_opt,liter,vcostl] = kappa_cost_optim(phifunct_n,g_nodal_n,theta_n,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n);
        [kappa_opt,liter,vcostl,remesh] = kappa_cost_only(phifunct_n,g_nodal_n,theta_n,theta_t,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,weights,h_C_0,Perim0,lambda_n,penalty_n,eta_n,iter,costfunc_n,vol_n,file_name,kappa_ini,Msmooth,Stiff_smooth,emass,flag_change_micro_ini,Group,file_gid,file_write,u0);
        
    case 'COST_ONLY_INC'
        [kappa_opt,liter,vcostl,remesh,iter_sub] = kappa_cost_only_remesh_iterative(phifunct_n,g_nodal_n,theta_n,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n,vol_n,file_name,kappa_ini,emass,iter_sub);
   
    case 'COST_BEST'
        [kappa_opt,liter,vcostl,remesh] =      kappa_best(phifunct_n,g_nodal_n,theta_n,theta_t,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,weights,h_C_0,Perim0,lambda_n,penalty_n,eta_n,iter,costfunc_n,vol_n,file_name,kappa_ini,Msmooth,Stiff_smooth,emass,flag_change_micro_ini,Group,file_gid,file_write,u0);
                                                   
end


end

