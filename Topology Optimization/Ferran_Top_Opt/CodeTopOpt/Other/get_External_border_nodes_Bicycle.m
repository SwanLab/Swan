function get_External_border_nodes_Bicycle (mesh,walls)
close all;clc;
if nargin == 0
    mesh = 'BicycleNew';
    walls = [1,2,3,4,5,6]; % 1 -> left // 2 -> right // 3 -> top // 4 -> bottom // 5 -> extra_hole // 6 -> Extra line between circle
end
coord = [];
eval(mesh);

% Plot all nodes
x = coord(:,2);
y = coord(:,3);
plot(x,y,'.k')
hold on

% Compute external border nodes
nnod = size(coord,1);
if any(walls == 1) % left wall in this case is a circle
    Center = [0.5,0.5];
    Radius = 0.5;
    eqn_circle = (x-Center(1)).^2 + (y-Center(2)).^2 - Radius^2;
    left_wall = (abs(eqn_circle) < 1e-3) & (x < Center(1));
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
if any(walls == 5)
    Center = [0.5,0.0999998];
    Radius = 0.65;
    eqn_circle = (x-Center(1)).^2 + (y-Center(2)).^2 - Radius^2;
    extra_hole1 = (abs(eqn_circle) < 1e-3) & (x > 0.8);
    
    Center = [0.5,0.0999998];
    Radius = 0.6;
    eqn_circle = (x-Center(1)).^2 + (y-Center(2)).^2 - Radius^2;
    extra_hole2 = (abs(eqn_circle) < 1e-3) & (x > 0.8);
    
    extra_hole = extra_hole1 | extra_hole2;
else
    extra_hole = false(nnod,1);
end
if any(walls == 6)
    m = 1.73073;
    n = -0.765352;
    line_vector = [1,m];
    point_line = [0,n];
    distance_point2line = @(x,y) (line_vector(1)*(y-point_line(2)) - line_vector(2)*(x-point_line(1)))/norm(line_vector);
    dist2line = distance_point2line(x,y);
    line = abs(dist2line) < 1e-3 & (x > 0.8) & (x < 0.83);
else
    line = false(nnod,1);
end


idx = left_wall | right_wall | top_wall | bottom_wall | extra_hole | line;
External_border_nodes = coord(idx,1);

% Plot external border nodes
x = coord(External_border_nodes,2);
y = coord(External_border_nodes,3);
plot(x,y,'xb')

% Generate code to add to the mesh
fprintf('ADD CODE IN CLIPBOARD TO MESH %sNew\n',mesh);
data = sprintf('External_border_nodes = [\n');
ext_nodes = sprintf('%i\n',External_border_nodes);
data = sprintf('%s%s];',data,ext_nodes);
clipboard('copy',data);
edit(mesh);

end