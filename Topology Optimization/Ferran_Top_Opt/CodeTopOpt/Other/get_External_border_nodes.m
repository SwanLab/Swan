function get_External_border_nodes (mesh,walls)
close all;clc;
if nargin == 0
    mesh = 'CantiliberbeamSymFineFineNew';
    walls = [1,2,3,4]; % 1 -> left // 2 -> right // 3 -> top // 4 -> bottom
end
gidcoord = [];
eval(mesh);

% Plot all nodes
x = gidcoord(:,2);
y = gidcoord(:,3);
plot(x,y,'.k')
hold on

% Compute external border nodes
nnod = size(gidcoord,1);
if any(walls == 1)
    left_wall = x == min(x);
else
    left_wall = false(nnod,1);
end
if any(walls == 2)
    right_wall = x == max(x);
else
    right_wall = false(nnod,1);
end
if any(walls == 3)
    top_wall = y == max(y);
else
    top_wall = false(nnod,1);
end
if any(walls == 4)
    bottom_wall = y == min(y);
else
    bottom_wall = false(nnod,1);
end
idx = left_wall | right_wall | top_wall | bottom_wall;
External_border_nodes = gidcoord(idx,1);

% Plot external border nodes
x = gidcoord(External_border_nodes,2);
y = gidcoord(External_border_nodes,3);
plot(x,y,'xb')

% Generate code to add to the mesh
fprintf('ADD CODE IN CLIPBOARD TO MESH %sNew\n',mesh);
data = sprintf('External_border_nodes = [\n');
ext_nodes = sprintf('%i\n',External_border_nodes);
data = sprintf('%s%s];',data,ext_nodes);
clipboard('copy',data);
edit(mesh);

end