function stop = outputfun_ipopt (iter,fval,data,outputopt)
stop = true;

optdata.fval = fval;
optdata.iteration = iter;
optdata.mu = data.mu;
optdata.incre_gamma = max(data.inf_pr,data.inf_du);
optdata.kappa = data.alpha_pr;
design_variable = data.x;
outputopt(design_variable,optdata,false);

% Use if IPOPT fails always at the same iteration. According to available
% information these crashes occur when IPOPT converges to machine precision
% https://bugs.launchpad.net/ubuntu/+source/coinor-ipopt/+bug/1388487

% if iter == 2605 
%     pause;
% end

end