function stop = plot_data_opt (design_variable,optdata,print_function,TYPE,verbosity,print_bool,plot_bool,last_iter,select_screen)
global post_info
global iter
stop = false; % used in FMINCON

%% Check data existance
% Volume
if isfield(post_info,'volume')
    volume = post_info.volume;
else
    volume = 0;
end
% Perimeter
if isfield(post_info,'perimeter')
    Perimeter = post_info.perimeter;
else
    Perimeter = 0;
end
% Compliance
if isfield(post_info,'compliance')
    compliance = post_info.compliance;
else
    compliance = 0;
end

% Cost
if isfield(optdata,'cost')
    cost = optdata.cost;
else
    if isfield(optdata,'fval') % check if using FMINCON
        cost = optdata.fval;
    else
        cost = 0;
    end
end
% Theta
if isfield(optdata,'theta')
    theta = optdata.theta*180/pi;
else
    if isfield(optdata,'mu')
        theta = optdata.mu;
    else
        theta = 0;
    end
end
% Kappa
if isfield(optdata,'kappa')
    kappa = optdata.kappa;
else
    if isfield(optdata,'stepsize') % check if using FMINCON
        kappa = optdata.stepsize;
    else
        kappa = 0;
    end
end
% Incre gamma
if isfield(optdata,'incre_gamma')
    incre_gamma = optdata.incre_gamma;
else
    incre_gamma = 0;
end
% Lambda
if isfield(optdata,'lambda')
    lambda = optdata.lambda;
    nconstr = length(lambda);
    if isfield(post_info,'lambda_fobj')
        if nconstr == 6
            lambda = insert_element(lambda,1,post_info.lambda_fobj);
        elseif nconstr == 12
            lambda = lambda(1:6); % MMA C-C* equality
            lambda = insert_element(lambda,1,post_info.lambda_fobj);
        else
            lambda = insert_element(lambda,2,post_info.lambda_fobj);
        end
    end
    nconstr = length(lambda);
    if nconstr < 2
        lambda(2) = 0;
    end
else
    lambda = [0,0];
end
% Penalty
if isfield(optdata,'penalty')
    penalty = optdata.penalty;
    nconstr = length(penalty);
    if nconstr < 2
        penalty(nconstr+1:2) = 0;
    end
else
    penalty = [0,0];
end
% Gradient
if isfield(optdata,'gradient')
    gradient = optdata.gradient;
else
    gradient = post_info.gradient;
end
% Last iter (correction for FMINCON)
if ~islogical(last_iter)
    last_iter = false;
end
% Iteration
if isfield(optdata,'iter_ls')
    iter_ls = optdata.iter_ls;
else
    if isfield(optdata,'iteration') % check if using FMINCON
        iter_ls = optdata.iteration;
    else
        iter_ls = 0;
    end
end

%% Post-process
% Save info and plot optimization data
if ~last_iter
    [post_info] = update_saved_variables(post_info,cost,theta,volume,lambda,Perimeter,kappa,incre_gamma,compliance,penalty(1),penalty(2),post_info.iepsilon,post_info.cumulative_iter + iter_ls);
end

% Print variables to GID
if print_bool || last_iter
    print_function(iter,post_info.gamma_reg_gp,post_info.gamma_nodal,gradient,design_variable,post_info.structural_values);
end

% Plots with optimization info
if (plot_bool && ~last_iter) || (last_iter && ~plot_bool)
    plot_variables(post_info.global_iter_n,post_info,last_iter);
end

% Verbosity
if ~last_iter && verbosity
    fprintf('\n');
    fprintf('incre_gamma = %g\n',incre_gamma);
    fprintf('theta = %g\n',theta);
    fprintf('volume = %g\n',volume);
    fprintf('perimeter = %g\n',Perimeter);
    fprintf('minmax = %g\t %g\n',min(post_info.gamma_reg_gp),max(post_info.gamma_reg_gp));
    fprintf('kappa = %g\n',kappa);
    fprintf('step = %g\n',post_info.iepsilon);
    
    if strcmp(TYPE,'MICRO')
        inv_matCh = inv(post_info.structural_values.matCh);
        nu12 = -inv_matCh(1,2)/inv_matCh(1,1);
        nu21 = -inv_matCh(2,1)/inv_matCh(2,2);
        E1 = 1/inv_matCh(1,1);
        E2 = 1/inv_matCh(2,2);
        fprintf('E1 = %16.16f\n',E1);
        fprintf('E2 = %16.16f\n',E2);
        fprintf('nu12 = %16.16f\n',nu12);
        fprintf('nu21 = %16.16f\n',nu21);
        fprintf('G = %16.16f\n',1/inv_matCh(3,3));
        matCh = post_info.structural_values.matCh;
        fprintf('Isotropic condition = %16.16f\n',matCh(1,1) - matCh(2,1) - 2*matCh(3,3));
        Ch_star_div = post_info.Ch_star;
        Ch_star_div (abs(Ch_star_div) < 1e-3) = 1;
        abs((post_info.structural_values.matCh - post_info.Ch_star)./Ch_star_div)
        if abs(nu12 - nu21) > 0.05
            pause;
        end
        
    end
    fprintf('\n');
end

if isempty(post_info.fhtri)
    fh = figure;
    mp = get(0, 'MonitorPositions');
    if size(mp,1) < select_screen
        select_screen = size(mp,1);
    end
    width = mp(1,3);
    height = mp(1,4);
    size_screen_offset = round([0.7*width,0.52*height,-0.71*width,-0.611*height],0);
    set(fh,'Position',mp(select_screen,:) + size_screen_offset);
    post_info.fhtri = trisurf(post_info.conectivities,post_info.coordinates(:,1),post_info.coordinates(:,2),double(post_info.gamma_nodal), ...
                              'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    view([0,90]);
    colormap(flipud(gray));
    set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
%     figure;
%     post_info.fhtri2 = trisurf(post_info.conectivities,post_info.coordinates(:,1),post_info.coordinates(:,2),gradient);
else
    set(post_info.fhtri,'FaceVertexCData',double(post_info.gamma_nodal));
    %colormap(flipud(gray))
    set(gca,'CLim',[0, 1],'XTick',[],'YTick',[]);
    
%     set(post_info.fhtri2,'FaceVertexCData',gradient);
    drawnow;
    %colorbar
end


end