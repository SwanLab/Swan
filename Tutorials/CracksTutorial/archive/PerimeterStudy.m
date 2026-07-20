%
clear;
clc;
close all;
%
%% Preprocess

% --- Intensity ---
fid = fopen('00001_intensity4.txt', 'r');
i = 1;
while ~feof(fid)
    if i<=7
        line = fgetl(fid);
        i = i+1;
    else
        line = fscanf(fid, '%lf');
        intensity = reshape(line,2,[])';
        break;
    end
end
fclose(fid);

% --- Nodes ---
fid = fopen('mesh_common5.txt', 'r');
i = 1;
while ~feof(fid)
    if i<=14
        line = fgetl(fid);
        i = i+1;
    else
        line = fscanf(fid, '%lf');
        coord = reshape(line,5,[])';
        break;
    end
end
fclose(fid);

% --- Elements ---
fid = fopen('mesh_common5.txt', 'r');
while ~feof(fid)
    line = fgetl(fid);
    if strcmp(strtrim(line), '[ELEMENTS]')
        fgetl(fid);   % skip the header
        break;
    end
end
line = fscanf(fid, '%lf');
connec = reshape(line, 5, [])';
fclose(fid);

coord     = coord(:,2:3);
connec    = connec(:,2:end);
intensity = intensity(:,2);

% Verification
% disp('Premier element :')
% disp(connec(1,:))
% disp('Coin 1 :'); disp(coord(connec(1,1),:))
% disp('Coin 2 :'); disp(coord(connec(1,2),:))
% disp('Coin 3 :'); disp(coord(connec(1,3),:))
% disp('Coin 4 :'); disp(coord(connec(1,4),:))

%% Mesh
s.coord  = coord;
s.connec = connec;
mesh = Mesh.create(s);

intFun = LagrangianFunction.create(mesh,1,'P1');
intFun.setFValues(intensity);
intFun.plot();

%% Filter
h = mesh.computeMeanCellSize();
e = 3*h;

sF.mesh       = mesh;
sF.filterType = 'PDE';
sF.trial      = intFun;
filter = Filter.create(sF);
filter.updateEpsilon(e);

intEpsFun = filter.compute(intFun,3);
intEpsFun.plot();

%% Perimeter computation
P = (2/e)*Integrator.compute(intFun.*(1-intEpsFun),mesh,2);