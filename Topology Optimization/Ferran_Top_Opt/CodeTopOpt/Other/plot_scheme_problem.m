function plot_scheme_problem (mesh)
close all;clc;
if nargin == 0
    mesh = 'TrencalosSupportFine';
end
coord = [];
External_border_nodes = [];
dirichlet_data = [];
pointload = [];
eval(mesh);
if exist('pointload_complete','var')
    pointload = pointload_complete;
end

% Plot all nodes
% x = coord(:,2);
% y = coord(:,3);
% plot(x,y,'.k')
% hold on

% Plot external border nodes
x = coord(External_border_nodes,2);
y = coord(External_border_nodes,3);
plot(x,y,'.b')
hold on

% Plot fixed nodes
x = coord(unique(dirichlet_data(:,1)),2);
y = coord(unique(dirichlet_data(:,1)),3);
plot(x,y,'+r')
hold on

% Plot loads
nodesolid = unique(pointload(:,1));
x = coord(nodesolid,2);
y = coord(nodesolid,3);
plot(x,y,'*g')
hold on


% Save image
fname = mesh;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
print(fname,'-dpng');

end