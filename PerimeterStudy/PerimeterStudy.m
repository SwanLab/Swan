
%
clear;
clc;
close all;
%

%% Preprocess

fid = fopen('00001_intensity4.txt', 'r');
i = 1;
while ~feof(fid)
    if i<=7
        line = fgetl(fid);
        disp(line)
        i = i+1;
    else
        line = fscanf(fid, '%lf');
        intensity = reshape(line,2,[])';
        break;
    end
end
fclose(fid);

fid = fopen('mesh_common4.txt', 'r');
i = 1;
while ~feof(fid)
    if i<=13
        line = fgetl(fid);
        disp(line)
        i = i+1;
    elseif i>14 && i<=51543
        line = fscanf(fid, '%lf');
        coord = reshape(line,5,[])';
        break;
    end
end
fclose(fid);

fid = fopen('mesh_common4.txt', 'r');
i = 1;
while ~feof(fid)
    if i<=103078
        line = fgetl(fid);
        disp(line)
        i = i+1;
    else
        line = fscanf(fid, '%lf');
        connec = reshape(line,5,[])';
        break;
    end
end
fclose(fid);

coord = coord(:,2:3);
connec = connec(:,2:end);
%connec = connec + 1;
intensity = intensity(:,2);

% coord(:,2) = -coord(:,2);
% coord(:,2) = coord(:,2) - min(coord(:,2)) + 0.5;

%% Mesh

s.coord = coord;
s.connec = connec;
mesh = Mesh.create(s);

intFun = LagrangianFunction.create(mesh,1,'P1');
intFun.setFValues(intensity);

intFun.plot();

%% Filter

h = mesh.computeMeanCellSize();
e = 3*h;

sF.mesh = mesh;
sF.filterType = 'PDE';
sF.trial = intFun;
filter = Filter.create(sF);
filter.updateEpsilon(e);

intEpsFun = filter.compute(intFun,3);

intEpsFun.plot();

%% Perimeter computation

P = (2/e)*Integrator.compute(intFun.*(1-intEpsFun),mesh,2);