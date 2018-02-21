function print_variables2gid_and_screen(iter,costfunc_n,theta_n,kappa_opt,vol_n,lambda,penalty,problembsc,d_u,file_gid,file_write,coordinates,element,...
        dim,gamma,g_nodal_n,gfunc_til_n,post,vdisp,nbdata,imesh,hC,constr,structural_values,Per,norm_g_n)

    switch problembsc.TYPE
    
        case 'MICRO'
            matCh = structural_values.matCh;
            
            fid = fopen([file_write,'Ch.txt'],'a+');
            switch problembsc.phisical_type
                case {'ELASTIC'}
                    fprintf(fid,'ITER: %3.0d  C11 %3.6f C12 %3.6f C13 %3.6f C22 %3.6f C23 %3.6f C33 %3.6f THETA  %10.6d VOL %10.6d hC %10.6d MESH %3.0f \n',...
                        iter,matCh(1,1),matCh(1,2),matCh(1,3),matCh(2,2),matCh(2,3),matCh(3,3),theta_n*180/pi,vol_n,hC,imesh);
                case {'THERMAL'}
                    fprintf(fid,'ITER: %3.0d  C11 %3.6f C12 %3.6f C22 %3.6f THETA  %10.6d VOL %10.6d hC %10.6d MESH %3.0f \n',...
                        iter,matCh(1,1),matCh(1,2),matCh(2,2),theta_n*180/pi,vol_n,hC,imesh);
            end
            
            
            
            fclose(fid);
            
    end
    
    fid = fopen([file_write,'executed_cases.txt'],'a+');
    fprintf(fid,'ITER: %3.1d  lambda %3.4f penalty %3.4f COST  %25.19d  THETA  %10.6d KAPPA %3.6f VOL %10.6d hC %10.6d Per %10.6d MESH %3.0f \n',...
        iter,lambda,penalty,costfunc_n,theta_n*180/pi,kappa_opt,vol_n,hC,Per,imesh);
    fclose(fid);
    printinfo(1,problembsc,d_u,file_gid,iter,coordinates,element,...
        dim,gamma,g_nodal_n,gfunc_til_n,structural_values,post,norm_g_n,vdisp,nbdata);
end
