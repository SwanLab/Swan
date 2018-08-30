function [gamma,lambda,fobj,iter] = projected_gradient(problem)

lambda = problem.lambda0;
penalty = problem.penalty;
gamma = problem.x0;
stopping_Criteria = 1;
iter = 1;
while stopping_Criteria
    
    gamma_old = gamma;
    stopping_criteria_line_search = 1;
    iter_ls = 1;
    
    printing = 1;
    gamma_print = gamma;
    [gamma_reg] = regularization_diff_reaction_equation(gamma,problem.epsilon,problem.Msmooth,problem.Stiff_smooth,problem.regularization);
    [gamma_gp,A_nod2gauss] = problem.nodal_2_elem(gamma_reg);
    [fobj(iter),gradient,~,vol(iter)] = lagrangian(problem,gamma_gp,lambda(iter),penalty,gamma_print,printing);
    gradient = regularization_diff_reaction_equation(A_nod2gauss'*gradient,problem.epsilon,problem.Msmooth,problem.Stiff_smooth,problem.regularization);
    
    t0 = norm(gamma_old)/norm(gradient);
    
     figure(1)
     hold on
     plot(problem.itert+iter,fobj(iter),'b+')
     drawnow
     
     figure(2)
     hold on
     plot(problem.itert+iter,vol(iter),'b+')
     drawnow
     
     figure(3)
     hold on
     plot(problem.itert+iter,lambda(iter),'b+')
     drawnow
     
   
     printing = 0;
    
    while stopping_criteria_line_search
        gamma_step = gamma_old - t0(iter_ls)*gradient;
        gamma = max(min(gamma_step,problem.ub),problem.lb);
        
        [gamma_reg] = regularization_diff_reaction_equation(gamma,problem.epsilon,problem.Msmooth,problem.Stiff_smooth,problem.regularization);
        
        [gamma_gp,A_nod2gauss] = problem.nodal_2_elem(gamma_reg);
        
        [fobj_new(iter_ls),gradient,constraint,volume(iter_ls)] = lagrangian(problem,gamma_gp,lambda(iter),penalty,gamma_print,printing);
        
        gradient = regularization_diff_reaction_equation(A_nod2gauss'*gradient,problem.epsilon,problem.Msmooth,problem.Stiff_smooth,problem.regularization);

        
        if fobj_new(iter_ls) < fobj(iter) || t0(iter_ls) < 1e-7
            stopping_criteria_line_search = 0;
        else 
            t0(iter_ls+1) = t0(iter_ls)/2;
        end
        
        
         
%         figure(4)
%         plot(t0(1:iter_ls),fobj_new(1:iter_ls),t0(1:iter_ls),fobj(iter)*ones(size(t0(1:iter_ls))))
%         drawnow
%         
        iter_ls = iter_ls + 1;
                
    end
    
    lambda(iter+1) = max(0,lambda(iter) + penalty*constraint);
    incre_gamma = (gamma - gamma_old)'*problem.Msmooth*(gamma - gamma_old)/(gamma_old'*problem.Msmooth*gamma_old)
    stopping_Criteria = incre_gamma >= 1e-6 || constraint > 1e-3;
    
    iter = iter + 1;
    

end
fobj_sol = fobj_new(end);
end

function [fobj,gradient,constraint,volume] = lagrangian(problem,gamma,lambda,penalty,gamma_print,printing)
[compliance,gradient_compliance] = problem.objective(gamma,gamma_print,printing);
constraint = problem.Aineq*gamma' - problem.bineq;
volume = problem.Aineq*gamma'*problem.Vfrac;
gradient_constraint = (problem.Aineq)';
    
fobj = compliance + lambda*constraint + 0.5*penalty*constraint^2;
gradient = gradient_compliance + lambda*gradient_constraint + penalty*constraint*gradient_constraint;
end