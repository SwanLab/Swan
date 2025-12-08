close all
clear all

% Specify the directory where the .mat files are located
directory = './AbrilTFGfiles/Data/20x20'; % Update this path as needed

% Get a list of all .mat files in the directory
files = dir(fullfile(directory, 'UL_*.mat'));

% Loop through each file and load it
i=1;
for k = 1:1:length(files)
    % Get the full path to the file
    filePath = fullfile(files(k).folder, files(k).name);

    % Load the file
    data=load(filePath);

    T(:,i) = data.T(:);  % This stores each file's contents in the cell array 'allData'
    Kcoarse(:,i) = data.K(:);

    disp(['Loaded: ', files(k).name]);  % Display the file being loaded
    i=i+1;
end


%xdata   = 1e-6:0.04:0.801;
xdata   = linspace(1e-6,0.999,20);
centers = xdata;
[fT,deim]   = parameterizedDataLagrange(T,xdata);
fK   = parameterizedDataLagrange(Kcoarse,xdata);

r=0.3;

u=fT(r);
z.mesh      = data.mesh;
z.order     = 'P1';
for i=1:8
  z.fValues   = reshape(u(:,i),[data.mesh.ndim,data.mesh.nnodes])';
  uFeFun = LagrangianFunction(z);%
  fileName = strrep("r" + num2str(r), '.', '_')+ "_SVDTraining" +num2str(i);
  centroids=computeCentroid(data.mesh);
  CoarsePlotSolution(uFeFun, data.mesh,[],fileName, r, centroids);
end



function [f,deim] = parameterizedDataLagrange(var,xdata)
deim    = DEIM(var);

coeff   = deim.rightVectors.';
s.coord = [xdata(1); xdata(end)];
s.connec = [1,2];
mesh = Mesh.create(s);
bFun = LagrangianFunction.create(mesh,1,'P8');

y =   leastSquares(coeff,xdata,bFun);
ss.mesh  = mesh;
ss.ndimf = 1;
ss.order = 'P8';

fR = arrayfun(@(i) ...
        LagrangianFunction(setfield(ss, 'fValues', y(:,i))), ...
        1:size(y,2), ...
        'UniformOutput', false);
fR = [fR{:}];

f = @(r) reshape( ...
        deim.basis * cell2mat(arrayfun(@(fun) evaluate(fun, 2*(r - min(xdata)) / (max(xdata) - min(xdata)) - 1), ...
        fR, 'UniformOutput', false)).', ...
        [], 8);

end


function y = leastSquares(data,xdata,f)
    xScaled = 2 * (xdata - min(xdata)) / (max(xdata) - min(xdata)) - 1;
    A = f.computeShapeFunctions(xScaled)';
    y = A \ data'; 
end

function centroids=computeCentroid(mesh)
    x0=mean(mesh.coord(:,1));
    y0=mean(mesh.coord(:,2));
    centroids = [x0,y0];
end