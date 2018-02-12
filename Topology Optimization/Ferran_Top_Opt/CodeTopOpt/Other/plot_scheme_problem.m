function plot_scheme_problem (mesh)
close all;clc;
if nargin == 0
    mesh = 'TrencalosSupportFine';
end
gidcoord = [];
External_border_nodes = [];
lnodes = [];
pointload = [];
eval(mesh);
if exist('pointload_complete','var')
    pointload = pointload_complete;
end

% Plot all nodes
% x = gidcoord(:,2);
% y = gidcoord(:,3);
% plot(x,y,'.k')
% hold on

% Plot external border nodes
x = gidcoord(External_border_nodes,2);
y = gidcoord(External_border_nodes,3);
plot(x,y,'.b')
hold on

% Plot fixed nodes
x = gidcoord(unique(lnodes(:,1)),2);
y = gidcoord(unique(lnodes(:,1)),3);
plot(x,y,'+r')
hold on

% Plot loads
nodesolid = unique(pointload(:,1));
x = gidcoord(nodesolid,2);
y = gidcoord(nodesolid,3);
plot(x,y,'*g')
hold on


% Save image
fname = mesh;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
print(fname,'-dpng');

end