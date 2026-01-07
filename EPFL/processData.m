close all
clear all
% Specify the directory where the .mat files are located
directory = './EPFL/dataLattice'; % Update this path as needed

% Get a list of all .mat files in the directory
files = dir(fullfile(directory, 'data_*.mat'));

% Loop through each file and load it
i=1;
for k = 1:4:length(files)
    % Get the full path to the file
    filePath = fullfile(files(k).folder, files(k).name);

    % Load the file
    load(filePath);

    % Optionally, store each loaded variable into a structure or cell array
    % For example, you could store the variables in a cell array:
    T(:,i) = EIFEoper.U;  % This stores each file's contents in the cell array 'allData'
    Td(:,i) = EIFEoper.Udef(:);
    Tr(:,i) = EIFEoper.Urb(:);
    Kcoase(:,i) = EIFEoper.Kcoarse(:);
    Kfine(:,i) = EIFEoper.Kfine(:);
    PhiD(:,i) = EIFEoper.PhiD(:);
    PhiR(:,i) = EIFEoper.PhiR(:);
    U(:,i) = EIFEoper.snapshots(:);
    % If you want to work with specific variables in each file, you can access them like:
    % data.variableName

    disp(['Loaded: ', files(k).name]);  % Display the file being loaded
    i=i+1;
end
xdata   = 0.01:0.04:0.961;
% N = 80;
% % Interval bounds
% a = 1e-6;
% b = 0.8;
% % Index vector
% i = 0:N;
% % Cosine spacing formula
% xdata = (a + b)/2 + (b - a)/2 * cos(pi * (1 - i / N));
centers = xdata;
[fT,deim,dfT]   = parameterizedDataLagrange(T,xdata);
% t    = fT(xdata);
[Tdef,~,dTdef]  = parameterizedDataLagrange(Td,xdata);
% def  = Tdef(xdata);
[Trb,~,dTrb]  = parameterizedDataLagrange(Tr,xdata);
% rb   = Trb(xdata);
[fK,~,dfK]   = parameterizedDataLagrange(Kcoase,xdata);
% Kcoarse = fK(xdata);

% for i=1:length(xdata)
%     error(i) = norm(def(:,i)+rb(:,i)-T(:,i))/norm(T(:,i));
%     error2(i) = norm(t(:,i)-T(:,i))/norm(T(:,i));

% end

% load('data_0.100.mat')
% EIFEoper.U = reshape(def,[],8) + reshape(rb',[],8);
% kcoarse = EIFEoper.U'*EIFEoper.Kfine*EIFEoper.U;
EIFEoper.Kcoarse = fK;
EIFEoper.Udef = Tdef;
EIFEoper.Urb  = Trb;
EIFEoper.U    = fT;
EIFEoper.dKcoarse = dfK;
EIFEoper.dUdef = dTdef;
EIFEoper.dUrb  = dTrb;
EIFEoper.dU    = dfT;

EIFEoper.deim    = deim;
filePath = './EPFL/parametrizedEIFEMLagrange20_der2_lattice.mat';
save(filePath,'EIFEoper')

% load('data_0.800.mat')
% def  = Tdef(0.1)lat
% rb   = Trb(0.1);
% % kcoarse = EIFEoper.U'*EIFEoper.Kfine*EIFEoper.U;
% % EIFEoper.Kcoarse = EIFEoper.U'*EIFEoper.Kfine*EIFEoper.U;;
% EIFEoper.Kcoarse = @(r) EIFEoper.Kcoarse;
% EIFEoper.Udef = @(r) def;
% EIFEoper.Urb  = @(r) rb;
% 
% filePath = '/home/raul/Documents/GitHub/EPFL/parametrizedEIFEM_T.mat';
% save(filePath,'EIFEoper')
% 
% EIFEoper.Kcoarse = @(r) (def)'*EIFEoper.Kfine*(def);
% EIFEoper.Udef = @(r) def;
% EIFEoper.Urb  = @(r) rb;
% 
% filePath = '/home/raul/Documents/GitHub/EPFL/parametrizedEIFEM_T_Kproj.mat';
% save(filePath,'EIFEoper')

load('./EPFL/data2/data_0.79723.mat')
% load('./EPFL/test/data_0.745.mat')
load('./EPFL/dataQ12/data_0.800_2.mat')
EIFEoper.Kcoarse = @(r) EIFEoper.Kcoarse;
EIFEoper.Udef = @(r) EIFEoper.Udef;
EIFEoper.Urb  = @(r) EIFEoper.Urb;

filePath = './EPFL/dataEIFEMQ12_2.mat';
save(filePath,'EIFEoper')

deim    = DEIM(var);
xdata   = 0.06:0.004:0.24;
centers = 0.06:0.004:0.24;
coeff   = deim.basis(deim.indices,:)\var(deim.indices,:);
f       = RBF(coeff',xdata,centers);

[Ut,St,Vt] = svd(T,"econ");
total = sum(diag(St).^2);
Splot = diag(St).^2/(total);
plot(Splot)

figure
[U,S,V] = svd(Kcoase,"econ");
total = sum(diag(S).^2);
Splot = diag(S).^2/(total);
plot(Splot)

function [f,deim] = parameterizedData(var,xdata,centers)
deim    = DEIM(var);

coeff   = deim.basis(deim.indices,:)\var(deim.indices,:);
% coeff = deim.rightVectors';
rbf       = RBF(coeff',xdata,centers);
% f = @ (r) deim.basis*rbf.evaluate(r) ;

% basis = deim.basis;
% for i = 1:size(deim.basis,2)
%     figure
%     plot(xdata,coeff(i,:))
%     hold on
%     plot(xdata,deim.rightVectors(:,i))
%     legend('MP','Singular vectors')
% end

f = @ (r) reshape( deim.basis*rbf.evaluate(r),[],8);

end

function [f,deim,df] = parameterizedDataLagrange(var,xdata)
deim    = DEIM(var);

coeff   = deim.basis(deim.indices,:)\var(deim.indices,:);
% dataFun = constructDataFEFunction2(coeff,xdata);
% m = dataFun(1).mesh;
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
% coeff = deim.rightVectors';

% for i=1:length(fR) 
%     yplot(:,i) = fR(i).evaluate(2*(xdata - min(xdata)) / (max(xdata) - min(xdata)) - 1); 
%     figure
%     plot(xdata,yplot(:,i))
%     hold on
%      plot(xdata,coeff(i,:))
% %     yplot2(i) = bFun.evaluate(xdata(i)); 
% end
% values = arrayfun(@(fun) evaluate(fun, 2*(0.1 - min(xdata)) / (max(xdata) - min(xdata)) - 1), fR, 'UniformOutput', false);
x_min = min(xdata);
L_dom = max(xdata) - x_min;

scale = 2 / L_dom;
shift = - (2 * x_min / L_dom) - 1;
B = deim.basis;
numFuns = length(fR);
f = @(r) reshape(B * evaluateAll(fR, r * scale + shift, numFuns), [], 8);

% f = @(r) reshape( ...
%         deim.basis * cell2mat(arrayfun(@(fun) evaluate(fun, 2*(r - x_min) / L_dom - 1), ...
%         fR, 'UniformOutput', false)).', ...
%         [], 8);
% 
gR_cell = arrayfun(@(fun) Grad(fun), fR, 'UniformOutput', false);
gR = [gR_cell{:}]; 
% --- Optimized Anonymous Function ---


% 3. Define the gradient function handle
% Note the multiplication by 'Jacobian' at the end of the inner calculation
% df = @(r) reshape( ...
%         deim.basis * (cell2mat(arrayfun(@(fun) evaluate(fun, 2*(r - x_min) / L_dom - 1), ...
%         gR, 'UniformOutput', false))'), ...
%         [], 8);

numGrads = length(gR);
df = @(r) reshape((B * evaluateAllGrads(gR, r * scale + shift, numGrads)), [], 8);

% f = @(r) deim.basis * cell2mat(arrayfun(@(fun) evaluate(fun, 2*(r - min(xdata)) / (max(xdata) - min(xdata)) - 1), ...
%         fR, 'UniformOutput', false)).';

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


% --- 2. Optimized Anonymous Function ---
% Note: We multiply the result by 'jac' to transform from local to global space


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

% function dataFun = constructDataFEFunction(coeff,xdata)
%     s.coord = xdata';
%     s.connec = [1:length(xdata)-1;2:length(xdata)]';
%     ss.mesh = Mesh.create(s);
%     ss.ndimf = 1;
%     ss.order = 'P1';
%     for i = 1:size(coeff,1)
%         ss.fValues = coeff(i,:)';
%         dataFun(i) = LagrangianFunction(ss);
%     end
% end

% function dataFun = constructDataFEFunction2(coeff, xdata)
%     % Build mesh once
%     s.coord  = xdata(:);                              % ensure row vector
%     s.connec = [(1:numel(xdata)-1)' (2:numel(xdata))']; % 2-column connectivity
%     ss.mesh  = Mesh.create(s);
%     ss.ndimf = 1;
%     ss.order = 'P1';
% 
%     % Construct one LagrangianFunction per row of coeff
%     dataFun = arrayfun(@(i) ...
%         LagrangianFunction(setfield(ss, 'fValues', coeff(i,:)')), ...
%         1:size(coeff,1), ...
%         'UniformOutput', false);
% 
%     % Convert cell array to vector of objects if class allows concatenation
%     if isa(dataFun{1}, 'LagrangianFunction')
%         dataFun = [dataFun{:}];
%     end
% end

function y = leastSquares(data,xdata,f)
    xScaled = 2 * (xdata - min(xdata)) / (max(xdata) - min(xdata)) - 1;
    A = f.computeShapeFunctions(xScaled)';
%     y =inv(A'*A)\A'*data';  
    y = A \ data'; 
end



% If you want to access specific variables in allData, use:
% allData{1}, allData{2}, ..., allData{k}
