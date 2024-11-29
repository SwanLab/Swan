function interpolateExample












end
% MATLAB Script to Read valores_interpolados.txt

% Specify the file name
filename = 'valores_interpolados.txt';

% Open the file and read the data
fileID = fopen(filename, 'r');

% Skip the header line (if it exists)
fgetl(fileID);  % Assuming the first line is a title

% Read the numerical data
data = fscanf(fileID, '%f');

% Close the file
fclose(fileID);

% Display the data
disp('Data from valores_interpolados.txt:');
disp(data)



% Specify the file name
filename = 'coordenadas_mesh2.txt';

fileContent = fileread(filename);

lines = regexp(fileContent, '\n', 'split');
numericData = strjoin(lines(2:end), '\n'); % Skip the header line

% Convert the numeric data into a matrix
data = sscanf(numericData, '%f,%f,%f');
data = reshape(data, 3, [])'; % Reshape into a Nx3 matrix
% Close the file
%fclose(fileID);



% Display the data
disp('Data from coordenadas_mesh2.txt:');
disp(data2)



