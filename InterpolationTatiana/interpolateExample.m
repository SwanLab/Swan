function interpolateExample
coord = readCoordinates();
alpha = alphaShape(coord);
plot(alpha);
[connec, vertices] = alphaTriangulation(alpha);
s.coord = coord;
s.connec = connec;
m = Mesh.create(s);

s.mesh = m;
s.order = 'P1';
s.fValues = readInterpolatedValues();
f = LagrangianFunction(s);
f.print('FirstData','Paraview')
end

function fValues = readInterpolatedValues()
    filename = 'valores_interpolados.txt';
    fileID = fopen(filename, 'r');
    fgetl(fileID); 
    fValues = fscanf(fileID, '%f');
    fclose(fileID);
end

function coord = readCoordinates()
    filename = 'coordenadas_mesh2.txt';
    fileContent = fileread(filename);
    lines = regexp(fileContent, '\n', 'split');
    numericData = strjoin(lines(2:end), '\n'); 
    data = sscanf(numericData, '%f,%f,%f');
    coord = reshape(data, 3, [])';
end


