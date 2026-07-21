clear;
clc;
close all;

%% Settings
MESH_FILE     = 'meshes/mesh_common.txt';
INTENSITY_DIR = 'meshes/intensities';
OUTPUT_CSV    = 'perimeters.csv';

%% Load shared geometry (once for all images)

% --- Nodes ---
fid = fopen(MESH_FILE, 'r');
i = 1;
while ~feof(fid)
    if i <= 14
        fgetl(fid);
        i = i + 1;
    else
        line = fscanf(fid, '%lf');
        coord = reshape(line, 5, [])';
        break;
    end
end
fclose(fid);

% --- Elements ---
fid = fopen(MESH_FILE, 'r');
while ~feof(fid)
    line = fgetl(fid);
    if strcmp(strtrim(line), '[ELEMENTS]')
        fgetl(fid);   % skip column header
        break;
    end
end
line   = fscanf(fid, '%lf');
connec = reshape(line, 5, [])';
fclose(fid);

coord  = coord(:, 2:3);
connec = connec(:, 2:end);

%% Build mesh and filter (once for all images)
s.coord  = coord;
s.connec = connec;
mesh = Mesh.create(s);

h = mesh.computeMeanCellSize();
e = 3 * h;

sF.mesh       = mesh;
sF.filterType = 'PDE';
% sF.trial      = LagrangianFunction.create(mesh, 1, 'P1');
filter        = Filter.create(sF);
filter.updateEpsilon(e);

%% Loop over positive and negative images
classes = {'positive', 'negative'};
labels  = [1, 0];

results = {};   % will hold {filename, label, perimeter}

for c = 1:2
    class_name  = classes{c};
    class_label = labels(c);
    class_dir   = fullfile(INTENSITY_DIR, class_name);

    files = dir(fullfile(class_dir, '*.txt'));

    fprintf('\n  Processing [%s] — %d files...\n', class_name, length(files));

    for k = 1:length(files)
        filepath = fullfile(class_dir, files(k).name);

        % --- Read intensity ---
        fid = fopen(filepath, 'r');
        i = 1;
        while ~feof(fid)
            if i <= 7
                fgetl(fid);
                i = i + 1;
            else
                line      = fscanf(fid, '%lf');
                intensity = reshape(line, 2, [])';
                break;
            end
        end
        fclose(fid);

        % --- Compute perimeter ---
        intFun = LagrangianFunction.create(mesh, 1, 'P1');
        intFun.setFValues(intensity(:, 2));
        intEpsFun = filter.compute(intFun, 3);
        P = (2/e) * Integrator.compute(intFun .* (1 - intEpsFun), mesh, 2);

        results{end+1} = {files(k).name, class_label, P};

        % Progress every 100 images
        if mod(k, 100) == 0 || k == length(files)
            fprintf('    %d/%d  last P = %.2f\n', k, length(files), P);
        end
    end
end

%% Save results to CSV
fid = fopen(OUTPUT_CSV, 'w');
fprintf(fid, 'filename,label,perimeter\n');
for i = 1:length(results)
    fprintf(fid, '%s,%d,%.6f\n', results{i}{1}, results{i}{2}, results{i}{3});
end
fclose(fid);

fprintf('\n  Saved %d results to %s\n', length(results), OUTPUT_CSV);