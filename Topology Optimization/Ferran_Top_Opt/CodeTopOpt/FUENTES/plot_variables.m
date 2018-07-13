function  plot_variables(iteration,post_info,last_iter)

global figure_handles
niter = length(iteration);
    
if last_iter || niter < 2 % initialize plots
    figure_handles = struct;
    figure_handles.cost = initialize_figure(1,'Cost','b');
    figure_handles.vol = initialize_figure(2,'Volume','b');
    figure_handles.lamV = initialize_figure(3,['\lambda_',strjoin(mat2cell(post_info.plot_names{1},1,ones(size(post_info.plot_names{1}))), '_')],'b');
    figure_handles.theta = initialize_figure(5,'\theta / \mu','b');
    figure_handles.kappa = initialize_figure(10,'\kappa','b','bar');
    figure_handles.per = initialize_figure(6,'Perimeter','b');
    figure_handles.incre_gamma = initialize_figure(9,'Incre gamma','b','semilogy');
    figure_handles.lamP = initialize_figure(7,['\lambda_',strjoin(mat2cell(post_info.plot_names{2},1,ones(size(post_info.plot_names{2}))), '_')],'b');
    figure_handles.comp = initialize_figure(11,'Compliance','b');
    figure_handles.ieps = initialize_figure(12,'nstep','b');
    figure_handles.rhoV = initialize_figure(4,['\rho_',strjoin(mat2cell(post_info.plot_names{1},1,ones(size(post_info.plot_names{1}))), '_')],'b');
    figure_handles.rhoP = initialize_figure(8,['\rho_',strjoin(mat2cell(post_info.plot_names{2},1,ones(size(post_info.plot_names{2}))), '_')],'b');
end   

% update plots
update_figure(figure_handles.cost,iteration,post_info.cost_n)
update_figure(figure_handles.vol,iteration,post_info.volume_n)
update_figure(figure_handles.theta,iteration,post_info.theta_n)
update_bar(figure_handles.kappa,iteration,post_info.kappa_n)
update_figure(figure_handles.per,iteration,post_info.Perimeter_n)
update_figure(figure_handles.incre_gamma,iteration,post_info.incre_gamma_n)
update_figure(figure_handles.comp,iteration,post_info.compliance_n)
update_figure(figure_handles.ieps,iteration,post_info.epsilon_n)
update_figure(figure_handles.rhoV,iteration,post_info.penalty_volume_n)
update_figure(figure_handles.rhoP,iteration,post_info.penalty_Perimeter_n)

[nconstr,~] = size(post_info.lambda_n);

% Lambda plots
if nconstr > 0
    if nconstr < 3 % st vol & per
        update_figure(figure_handles.lamV,iteration,post_info.lambda_n(1,:))
        update_figure(figure_handles.lamP,iteration,post_info.lambda_n(2,:))
    elseif nconstr >= 3 && nconstr < 7 % st enforce_Ch_diff
        % Creation of new line handles
        if length(figure_handles.lamV) < nconstr
            figure_handles.axh_constr = figure_handles.lamV.Parent;
            figure_handles.lamV = {figure_handles.lamV};
            hold on
            line_colour = {'b','r','g','y','k','m'};
            for i = 2:nconstr
                figure_handles.lamV{i} = add_lines_plot (figure_handles.axh_constr,line_colour{i});
            end
        end
        
        % Update plot
        for i = 1:nconstr
            update_figure(figure_handles.lamV{i},iteration,post_info.lambda_n(i,:))
        end
        legend(figure_handles.axh_constr,post_info.constr_names,'Location','northwest');
    else % st per, vol & enforce_Ch_diff
        idx = nconstr-6; % show 2 constraints in 1st plot and the rest in 2nd
        % Creation of new line handles
        first_time = (length(figure_handles.lamP) + length(figure_handles.lamV)) < nconstr;
        if first_time
            figure_handles.axh_constr2 = figure_handles.lamP.Parent;
            figure_handles.lamP = {figure_handles.lamP};
            hold on
            line_colour = {'b','r','g','y','k','m'};
            for i = 2:nconstr-idx
                figure_handles.lamP{i} = add_lines_plot (figure_handles.axh_constr2,line_colour{i});
            end
            figure_handles.axh_constr1 = figure_handles.lamV.Parent;
            figure_handles.lamV = {figure_handles.lamV};
            hold on
            line_colour = {'b','r','g','y','k','m'};
            for i = 1:idx
                figure_handles.lamV{i} = add_lines_plot (figure_handles.axh_constr1,line_colour{i});
            end
        end
        
        % Update plot
        for i = 1:idx
            update_figure(figure_handles.lamV{i},iteration,post_info.lambda_n(i,:))
        end
        for i = (idx+1):nconstr
            update_figure(figure_handles.lamP{i-idx},iteration,post_info.lambda_n(i,:))
        end
        legend(figure_handles.axh_constr1,post_info.constr_names1);
        legend(figure_handles.axh_constr2,post_info.constr_names2,'Location','northwest');
    end
end

end

function handle = initialize_figure(number,tit,style,type)

if nargin < 4
    type = 'plot';
end

subplot(3,4,number)
switch type
    case 'bar'
        handle = bar(0,0,style);
    case 'semilogy'
        handle = semilogy(0,0,style);
    case 'plot'
        handle = plot(0,0,style);
        
end
title(tit);

end

function update_figure(h,iteration,variable)

    set(h,'XData',iteration,'YData',variable);
    drawnow

end

function update_bar(h,iteration,variable)

    set(h,'XData',iteration,'YData',variable);
    drawnow

end

function fh_new = add_lines_plot (fh_parent,lcolor)
axes(fh_parent);
fh_new = line;
% fh_new = line(fh_parent);
set(fh_new,'Color',lcolor);

end