function [kappa_opt,liter,vcostl,remesh,iter_sub] = kappa_cost_only_remesh(phifunct_l,g_nodal_n,theta_n,norm_g_n,element,dim,coordinates,problembsc,fixnodes,fext,iter,costfunc_n,vol_n,file_name,kappa_ini,emass,iter_sub)

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

    while (liter <= lmax && stop==0)
    
               
        [phifunct_l1] = update_phifunc(theta_n,kappa_l,phifunct_l,g_nodal_n,norm_g_n);
        [phigp_l1] = interpol(phifunct_l1,element,dim,problembsc);
        
        [vol_l1] = cal_omega(phifunct_l1,dim,element,problembsc,coordinates);
        
        constr = call_constraint(vol_l1,element,problembsc);
        % solve ku=f
        [~,matCh] = module_M(phifunct_l1,phigp_l1,element,fixnodes,problembsc,coordinates,fext,dim);
        
        %element1 = update_lambda(problembsc,element,vol_l1);
        %[cost_l1] = cal_cost_funct(matCh,problembsc,element1,fext,vol_l1);
        
        % cost function
        [cost_l1] = cal_cost_funct(matCh,problembsc,element,fext,vol_l1);
        
         
        [theta_l1] = cal_theta(dim,element,g_nodal_n,phifunct_l1,emass);
         
        
        
        vcostl(liter) = cost_l1;     
        kappa(liter) = kappa_l;
  
%         figure(11)
%         plot(kappa(1:liter),vcostl(1:liter),'-+')
%         hold on
%         plot([0 1],[costfunc_n costfunc_n],'r')
%         hold off

           % fid = fopen([file_name,'executed_cases.txt'],'a+');
   
           % fprintf(fid,'ITER: %3.0d LITER: %3.0d COST_L %25.19d THETA  %10.6d KAPPA %3.6f VOL %10.6d \n',...
           %     iter,liter,cost_l1,theta_n*180/pi,kappa_l,vol_l1);
           % fclose(fid);
   

           
            
        if (cost_l1<costfunc_n*(1+extra_tol)) && abs(vol_l1-vol_n)/vol_n <= element.Vol_inc_tol;%0.05;%*(theta_n/(2*pi/3))^2)) %|| (theta_l1 < theta_n)
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
    if (liter > lmax)  % maximo numero de iteraciones
        
        switch problembsc.lagrange_update;
            case 'LINEAL'
                criterio = element.imesh < 7;
                
            case 'AUGMENTED'
                criterio = element.imesh < 7 && (abs(constr) <= 8*1e-2);
                
        end
   
        
        if iter_sub == 1
            remesh = 1;
            kappa_opt = kappa_min;
            iter_sub = 0;
        else
            remesh = 0;
            kappa_opt = 0;
            iter_sub = iter_sub + 1;
        end
        
                
                
    else
        iter_sub = 0;
        remesh = 0;
%         if criterio
%         remesh = 1;
%         else
%         remesh = 0;    
%         end
       
        %kappa_opt = kappa_min;

    end
    
end



    

