function plot_check_derivative(g,g0,gv,gv0,gp,gp0,coordinates,Msmooth,save_fig)

%% Field plots
fname = 'Compliance Field';
tit{1} = 'Finite Differences';
tit{2} = 'User-supplied';
plot_field(fname,g,g0,coordinates,tit,save_fig);

fname = 'Volume Field';
tit{1} = 'Finite Differences';
tit{2} = 'User-supplied';
plot_field(fname,gv,gv0,coordinates,tit,save_fig);

fname = 'Perimeter Field';
tit{1} = 'Finite Differences';
tit{2} = 'User-supplied';
plot_field(fname,gp,gp0,coordinates,tit,save_fig);

%% Gradient differences
fname = 'Compliance Difference';
tit = 'Finite differences minus user-supplied';
plot_figure(fname,g,g0,coordinates,tit,Msmooth,save_fig);

fname = 'Volume Difference';
tit = 'Finite differences minus user-supplied';
plot_figure(fname,gv,gv0,coordinates,tit,Msmooth,save_fig);

fname = 'Perimeter Difference';
tit = 'Finite differences minus user-supplied';
plot_figure(fname,gp,gp0,coordinates,tit,Msmooth,save_fig);



end

function plot_field(fname,g,g0,coordinates,tit,save_fig)

fh = figure('Name',fname);
subplot(1,2,1);
plot_nodal_field(g,coordinates);
title(tit{1});
subplot(1,2,2);
plot_nodal_field(g0,coordinates);
title(tit{2});
if save_fig
    savefig(fh,fname);
    print(fh,fname,'-dpng');
    close(fh);
end

end

function plot_figure(fname,g,g0,coordinates,tit,Msmooth,save_fig)

fh = figure('Name',fname);
% plot(g-g0,'x');
een = ones(size(g));
factor = (g'*Msmooth*g)/(een'*Msmooth*een);
plot_nodal_field((g-g0)/factor,coordinates);
title(tit);
if save_fig
    savefig(fh,fname);
    print(fh,fname,'-dpng');
    close(fh);
end

end

