function [post_info] = update_saved_variables(post_info,cost,theta,volume,lambda,Perimeter,kappa,incre_gamma,compliance,penalty_volume,penalty_Perimeter,epsilon,global_iter)   
post_info.cost_n = [post_info.cost_n cost];
post_info.theta_n = [post_info.theta_n theta];
post_info.volume_n = [post_info.volume_n volume];
post_info.lambda_n = [post_info.lambda_n lambda(:)];
post_info.Perimeter_n = [post_info.Perimeter_n Perimeter];
post_info.kappa_n = [post_info.kappa_n kappa];
post_info.incre_gamma_n = [post_info.incre_gamma_n incre_gamma];
post_info.compliance_n = [post_info.compliance_n compliance];
post_info.penalty_volume_n = [post_info.penalty_volume_n penalty_volume];
post_info.penalty_Perimeter_n = [post_info.penalty_Perimeter_n penalty_Perimeter];
post_info.epsilon_n = [post_info.epsilon_n epsilon];
post_info.global_iter_n = [post_info.global_iter_n global_iter];
end