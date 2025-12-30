close all
clear all

% Specify case parameters
p.nelem     =  20;
p.Inclusion = 'Material';         % 'Hole'/'Material'/'HoleRaul'
p.Sampling  = 'Isolated';     % 'Isolated'/'Oversampling'
meshName    = p.nelem+"x"+p.nelem;

% Specify the directory where the .mat files are located

% Get a list of all .mat files in the directory
files = dir(fullfile("AbrilTFGfiles/Data/",p.Inclusion,p.Sampling,meshName, 'r0_*.mat'));

% Loop through each file and load it
i=1;
for k = 1:1:length(files)
    % Get the full path to the file
    filePath = fullfile(files(k).folder, files(k).name);

    % Load the file
    load(filePath,'EIFEoper',"mesh");

    T(:,i) = EIFEoper.U(:);  % This stores each file's contents in the cell array 'allData'
    Td(:,i) = EIFEoper.Udef(:);
    Tr(:,i) = EIFEoper.Urb(:);
    Kcoarse(:,i) = EIFEoper.Kcoarse(:);
    Kfine(:,i) = EIFEoper.Kfine(:);
    PhiD(:,i) = EIFEoper.PhiD(:);
    PhiR(:,i) = EIFEoper.PhiR(:);

    disp(['Loaded: ', files(k).name]);  % Display the file being loaded
    i=i+1;
end

%xdata   = linspace(0,0.9,10);
xdata   = linspace(0,0.95,20);
%xdata   = linspace(0,0.98,50);

centers = xdata;
%mesh=data.mesh;

[fT,deim,dfT,fR]   = parameterizedDataLagrange(T,xdata);
[Tdef,~,dTdef]  = parameterizedDataLagrange(Td,xdata);
[Trb,~,dTrb]  = parameterizedDataLagrange(Tr,xdata);
[fK,~,dfK,~]   = parameterizedDataLagrange(Kcoarse,xdata);

fileName=fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,meshName,"HOfunction.mat");
save(fileName,"fR","fT","deim","mesh");

EIFEoper.Kcoarse    = fK;
EIFEoper.Udef       = Tdef;
EIFEoper.Urb        = Trb;
EIFEoper.U          = fT;
EIFEoper.dKcoarse   = dfK;
EIFEoper.dUdef      = dTdef;
EIFEoper.dUrb       = dTrb;
EIFEoper.dU         = dfT;
EIFEoper.deim       = deim;


filePath = fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,meshName,"parametrizedEIFEM.mat");
save(filePath,'EIFEoper');


function [f,deim,df,fR] = parameterizedDataLagrange(var,xdata)
deim    = DEIM(var);
coeff   = deim.basis(deim.indices,:)\var(deim.indices,:);
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

x_min = min(xdata);
L_dom = max(xdata) - x_min;

scale = 2 / L_dom;
shift = - (2 * x_min / L_dom) - 1;
B = deim.basis;
numFuns = length(fR);

f = @(r) reshape(B * evaluateAll(fR, r * scale + shift, numFuns), [], 8);

gR_cell = arrayfun(@(fun) Grad(fun), fR, 'UniformOutput', false);
gR = [gR_cell{:}]; 

numGrads = length(gR);
df = @(r) reshape((B * evaluateAllGrads(gR, r * scale + shift, numGrads)), [], 8);

end


% --- Helper Function (Place at end of script or in separate file) ---
function vals = evaluateAll(fR, s, n)
    % Pre-allocate the matrix: rows = functions, cols = points
    % This assumes evaluate returns a scalar or row vector per point
    numPoints = length(s);
    vals = zeros(n, numPoints);
    for i = 1:n
        vals(i, :) = evaluate(fR(i), s);
    end
end



% --- Helper Function (Place at end of script) ---
function vals = evaluateAllGrads(gR, s, n)
    % Pre-allocate matrix: [Number of Basis Functions x Number of Points]
    numPoints = length(s);
    vals = zeros(n, numPoints);
    for i = 1:n
        % Direct assignment is significantly faster than cell2mat
        vals(i, :) = evaluate(gR(i), s);
    end
end


function y = leastSquares(data,xdata,f)
    xScaled = 2 * (xdata - min(xdata)) / (max(xdata) - min(xdata)) - 1;
    A = f.computeShapeFunctions(xScaled)';
    y = A \ data'; 
end

