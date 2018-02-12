function [kappa_opt,liter,vcostl,remesh,iter_sub] = kappa_cost_only_remesh_iterative(phifunct_l,g_nodal_n,theta_n,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n,vol_n,file_name,kappa_ini,emass,iter_sub)

freduc = element.material.kappa_reduction;

kappa_min = element.material.kappa_min;
kappa_end = element.material.kappa_end;

lmax = element.material.kappa_maxiter;
liter = 1; 
switch kappa_ini
    case 0
        kappa_l = 1;    
        extra_tol = element.fobj_inc_tol*((theta_n/(2*pi/3))^1);
    case kappa_ini
        kappa_l = 1;    
        extra_tol = 0;
        %extra_tol = 1e2*((theta_n/(2*pi/3))^1);
    otherwise
        kappa_l = min(1,1.5*kappa_ini);
        %kappa_l = 1;
        extra_tol = 0;
        %extra_tol = 1e2*((theta_n/(2*pi/3))^1);
end

stop = 0;

%fid = fopen([file_name,'executed_cases.txt'],'a+');
%fprintf(fid,'ITER: %3.0d            COST_L %25.19d THETA  %10.6d  \n',iter,costfunc_n,theta_n*180/pi);
%fclose(fid);
cost_l1 = 1+costfunc_n;
%kappa = 1./(2.^[0:10])';
kappa = 1./(2.^(0.5*[0:20]))';

lmax = length(kappa);
vcostl = zeros(lmax,1);

    while cost_l1 > costfunc_n && liter <= lmax
         kappa_l = kappa(liter);
         cost_l1 = call_new_cost(kappa_l,theta_n,phifunct_l,g_nodal_n,norm_g_n,element,dim,problembsc,fixnodes,coordinates,fext);
         vcostl(liter,1) = cost_l1;     
         kappa_opt = kappa_l;
         liter = liter +1;
         remesh = 0;
    end
    
    if cost_l1 > costfunc_n
    liter = 0;
    inc_fun_obj = 0.01;
    cost = 1+costfunc_n*(1+inc_fun_obj);
    
    while cost > costfunc_n*(1+inc_fun_obj)
        liter = liter + 1; 
        cost = vcostl(liter);
       
       if liter >= lmax
           liter = 0;
           inc_fun_obj = inc_fun_obj*2;
           cost = 1+costfunc_n*(1+inc_fun_obj);
       end
    end
    inc_fun_obj
    kappa_opt = kappa(liter)
    end
end



function cost_l1 = call_new_cost(kappa_l,theta_n,phifunct_l,g_nodal_n,norm_g_n,element,dim,problembsc,fixnodes,coordinates,fext)
        [phifunct_l1] = update_phifunc(theta_n,kappa_l,phifunct_l,g_nodal_n,norm_g_n);
        [phigp_l1] = interpol(phifunct_l1,element,dim,problembsc);
        
        [vol_l1] = cal_omega(phifunct_l1,dim,element,problembsc,coordinates);
        
        %constr = call_constraint(vol_l1,element,problembsc);
        % solve ku=f
        [~,matCh] = module_M(phifunct_l1,phigp_l1,element,fixnodes,problembsc,coordinates,fext,dim);
        
        %element1 = update_lambda(problembsc,element,vol_l1);
        %[cost_l1] = cal_cost_funct(matCh,problembsc,element1,fext,vol_l1);
        
        % cost function
        [cost_l1] = cal_cost_funct(matCh,problembsc,element,fext,vol_l1);
end

