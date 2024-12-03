function interpolateExample
[coord,fValues] = readData();%alpha = alphaShape(coord);
%plot(alpha);
%[connec] = delaunayTriangulation(coord);
points = coord;

%rng(1); % Seed for reproducibility
%numPoints = 200;
%points = rand(numPoints, 3); % Random 3D points



%boundaryIndices = boundary(points(:,1), points(:,2), 0.8); % Adjust shrink factor (0-1)



% Create an alpha shape with a larger initial alpha value
alpha = 1.0; % Start with a larger alpha
shp = alphaShape(points, alpha);

% Check if the alpha shape is empty
while isempty(shp.Points)
    alpha = alpha * 1.5; % Increase alpha and try again
    shp = alphaShape(points, alpha);
end

% Delaunay Triangulation in 3D
dt = delaunayTriangulation(points);

% Get the tetrahedra (connectivity list)
tetrahedra = dt.ConnectivityList;

% Compute the centroids of each tetrahedron
centroids = (points(tetrahedra(:,1), :) + ...
             points(tetrahedra(:,2), :) + ...
             points(tetrahedra(:,3), :) + ...
             points(tetrahedra(:,4), :)) / 4;

% Check which centroids are inside the alpha shape
inside = inShape(shp, centroids);

% Keep only tetrahedra whose centroids are inside the alpha shape
filteredTetrahedra = tetrahedra(inside, :);

% Extract unique triangular faces from filtered tetrahedra
faces = freeBoundary(triangulation(filteredTetrahedra, points));

% Plot the results
figure;
trisurf(faces, points(:,1), points(:,2), points(:,3), ...
        'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Non-convex boundary mesh
hold on;
plot3(points(:,1), points(:,2), points(:,3), 'k.', 'MarkerSize', 10); % Scattered points
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Non-Convex 3D Mesh from Scattered Points');
view(3);


s.coord = coord;
s.connec = filteredTetrahedra;
m = Mesh.create(s);

s.mesh = m;
s.order = 'P1';
s.fValues = fValues;
f = LagrangianFunction(s);
f.print('FirstData','Paraview')
end

function [coord,fValues] = readData()
    filename = 'FDJ.txt';
    fileContent = fileread(filename);
    lines = regexp(fileContent, '\n', 'split');
    numericData = strjoin(lines(2:end), '\n'); 
    data = sscanf(numericData, '%f,%f,%f,%f');
    dataV = reshape(data, 4, [])';
    coord = dataV(:,1:3);
    fValues = dataV(:,4);
end


