function [stopping_criteria_opt,kernel_case] = compute_stopping_criteria(theta,incre_gamma,constraint,problem,kernel_case,algorithm,iepsilon, penalty_volume_0)
global iter
switch algorithm
    case 'level_set'
        
        if penalty_volume_0 == 0
            stopping_criteria_opt = theta >= problem.theta_min(iepsilon) || iter < 10;
            %change_kernel = theta <= problem.theta_min_P1 || abs(constraint) <= problem.constraint_tol(iepsilon);
        else
            stopping_criteria_opt = theta >= problem.theta_min(iepsilon) || abs(constraint) >= problem.constraint_tol(iepsilon) || iter < 10;    
        end
        
        %stopping_criteria_opt = theta >= problem.theta_min(iepsilon) || iter < 10;
        %stopping_criteria_opt = incre_gamma >= problem.incre_gamma_min || iter < 10;
        %stopping_criteria_opt = incre_gamma >= problem.incre_gamma_min || abs(constraint) >= problem.constraint_tol;
        
        
   %     if change_kernel && strcmp(kernel_case,'P1_kernel')
   %         stopping_criteria_opt = 1;
   %         kernel_case = 'P0_kernel';
   %     end
        
    case 'Projected_gradient'
        
        if penalty_volume_0 == 0
            stopping_criteria_opt = incre_gamma >= problem.incre_gamma_min(iepsilon) || iter < 10;
            %change_kernel = theta <= problem.theta_min_P1 || abs(constraint) <= problem.constraint_tol(iepsilon);
        else
            stopping_criteria_opt = incre_gamma >= problem.incre_gamma_min(iepsilon) || abs(constraint) >= problem.constraint_tol(iepsilon) || iter < 10;
        end
        
        %stopping_criteria_opt = incre_gamma >= problem.incre_gamma_min || abs(constraint) <= problem.constraint_tol(iepsilon);
%        stopping_criteria_opt = incre_gamma >= problem.incre_gamma_min || abs(constraint) <= problem.constraint_tol(iepsilon);
        %stopping_criteria_opt = incre_gamma >= problem.incre_gamma_min || abs(constraint) >= problem.constraint_tol;
       % change_kernel = incre_gamma <= problem.incre_gamma_min_P1 || abs(constraint) <= problem.constraint_tol;
        
   %     if change_kernel && strcmp(kernel_case,'P1_kernel')
   %         stopping_criteria_opt = 1;
   %         kernel_case = 'P0_kernel';
   %     end
        
        
end

end