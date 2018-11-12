function target_parameters = update_target_parameters(target_parameters,t)

target_parameters.Cstar = update_parameter(target_parameters.Cstar_ini,target_parameters.Cstar_final,target_parameters.alpha_Cstar(t));
target_parameters.Vtarget = update_parameter(target_parameters.vol_ini,target_parameters.vol_final,target_parameters.alpha_vol(t));
target_parameters.epsilon = update_parameter(target_parameters.epsilon_ini,target_parameters.epsilon_final,target_parameters.alpha_eps(t));
target_parameters.epsilon_isotropy = update_parameter(target_parameters.epsilon_isotropy_ini,target_parameters.epsilon_isotropy_final,target_parameters.alpha_isotropy2d(t));

target_parameters.PrimalTol = update_parameter(target_parameters.PrimalTol_ini,target_parameters.PrimalTol_final,target_parameters.alpha_PrimalTol(t));
target_parameters.DualTol = update_parameter(target_parameters.DualTol_ini,target_parameters.DualTol_final,target_parameters.alpha_DualTol(t));

end

function x = update_parameter(x_ini,x_final,alpha)

x = (1-alpha)*x_ini + alpha*x_final;

end