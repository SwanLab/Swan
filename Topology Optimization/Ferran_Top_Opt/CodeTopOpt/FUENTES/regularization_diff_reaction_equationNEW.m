function [f_reg,f_nodal,f_gp] = regularization_diff_reaction_equation(f,epsilon,kernel_case,Msmooth,Stiff_smooth,regularization,coordinates,element,problembsc,dim,P_operator,A_nodal_2_gauss,variable_case,M0,smoothing_function,algorithm)
%Msmooth = diag(sum(Msmooth));


switch kernel_case
    case 'PDE'
        switch variable_case
            case 'gamma'
                switch  algorithm
                    
                    case 'level_set'
                        M2 = faireF2(coordinates',element.conectivities',f);
                        Rinv  = (epsilon^2*Stiff_smooth + Msmooth);
                        dim.nunkn = 1;
                        f_reg = solverp(zeros(size(M2)),Rinv,element.pnods,M2,[],dim);%f_reg = Rinv\M2;
                        f_nodal = f < 0;
                        f_gp = A_nodal_2_gauss*f_reg;

                        
                        
                    case 'Projected_gradient'
                        rhs = Msmooth*f;
                        Rinv  = (epsilon^2*Stiff_smooth + Msmooth);
                        dim.nunkn = 1;
                        f_reg = solverp(zeros(size(rhs)),Rinv,element.pnods,rhs,[],dim);%f_reg = Rinv\rhs;
                        f_nodal = f;
                        f_gp = A_nodal_2_gauss*f_reg;
                        
                end
                
            case 'gradient'
                switch  algorithm
                    
                    case 'level_set'
                        rhs = (A_nodal_2_gauss'*M0*f);
                        Rinv = (epsilon^2*Stiff_smooth + Msmooth);
                        dim.nunkn = 1;
                        f_reg = solverp(zeros(size(rhs)),Rinv,element.pnods,rhs,[],dim);%f_reg = Rinv\rhs;
                        f_nodal = f_reg;
                        f_gp = A_nodal_2_gauss*f_reg;
                        
                    case 'Projected_gradient'
                        rhs = (A_nodal_2_gauss'*M0*f);
                        Rinv = (epsilon^2*Stiff_smooth + Msmooth);
                        dim.nunkn = 1;
                        f_reg = solverp(zeros(size(rhs)),Rinv,element.pnods,rhs,[],dim);%f_reg = Rinv\rhs;
                        
%                         rhs = (A_nodal_2_gauss'*f);
%                         Rinv = (epsilon^2*Stiff_smooth + Msmooth);
%                         f_reg = Msmooth*(Rinv\rhs);
                        
                        
                        
                        f_nodal = f_reg;
                        f_gp = A_nodal_2_gauss*f_reg;
                end
                
                
        end
        
        
    case 'P1_kernel'
        
        switch variable_case
            case 'gamma'
                switch  algorithm
                    
                    case 'level_set'
                        M2 = faireF2(coordinates',element.conectivities',f);
                        f_gp = P_operator*M2;
                        f_reg = A_nodal_2_gauss'*f_gp;
                        f_nodal = f < 0;
                        
                    case 'Projected_gradient'
                        rhs = Msmooth*f;
                        f_gp = P_operator*rhs;
                        f_reg = f;
                        f_nodal = f_reg;
                        
                        
                end
                
            case 'gradient'
                switch  algorithm
                    
                    case 'level_set'
                        f_reg = P_operator'*M0*f;
                        f_nodal = 0;
                        f_gp = 0;
                        
                    case 'Projected_gradient'
                        %f_reg = P_operator'*M0*f;
                        %f_reg = Msmooth*(P_operator'*f);
                        f_reg = P_operator'*M0*f;
                        f_nodal = 0;
                        f_gp = 0;
                end
                
                
        end
        
        
    case 'P0_kernel'
        
        switch variable_case
            case 'gamma'
                switch  algorithm
                    
                    case 'level_set'
                        f_gp = cal_vol_mat(f,dim,element,problembsc,coordinates)';
                        
                        %M2 = faireF2(coordinates',element.conectivities',f);
                        %f_gp = A_nodal_2_gauss*(Msmooth\(M2));
                        
                        f_reg = A_nodal_2_gauss'*f_gp;
                        f_nodal = f < 0;
                        
                    case 'Projected_gradient'
                        f_gp = A_nodal_2_gauss*f;
                        f_reg = f;
                        f_nodal = f;
                end
                
            case 'gradient'
                switch  algorithm
                    
                    case 'level_set'
                        %f_reg = A_nodal_2_gauss'*f;
                        f_reg = Msmooth\(A_nodal_2_gauss'*M0*f); 
                        f_nodal = 0;
                        f_gp = 0;
                        
                    case 'Projected_gradient'
                        f_reg = Msmooth\(A_nodal_2_gauss'*M0*f); 
                       % f_reg = (A_nodal_2_gauss'*f); 
                        f_nodal = 0;
                        f_gp = 0;
                end
                
                
        end
        
        
end
        

end