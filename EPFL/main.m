close all
clear all
% Specify the directory where the .mat files are located
directory = './EPFL/'; % Update this path as needed

% Get a list of all .mat files in the directory
files = dir(fullfile(directory, 'data_*.mat'));

% Loop through each file and load it
for k = 1:length(files)
    % Get the full path to the file
    filePath = fullfile(files(k).folder, files(k).name);

    % Load the file
    load(filePath);

    % Optionally, store each loaded variable into a structure or cell array
    % For example, you could store the variables in a cell array:
    T(:,k) = EIFEoper.U;  % This stores each file's contents in the cell array 'allData'
    Td(:,k) = EIFEoper.Udef(:);
    Tr(:,k) = EIFEoper.Urb(:);
    Kcoase(:,k) = EIFEoper.Kcoarse(:);
    Kfine(:,k) = EIFEoper.Kfine(:);
    PhiD(:,k) = EIFEoper.PhiD(:);
    PhiR(:,k) = EIFEoper.PhiR(:);
    U(:,k) = EIFEoper.snapshots(:);
    % If you want to work with specific variables in each file, you can access them like:
    % data.variableName

    disp(['Loaded: ', files(k).name]);  % Display the file being loaded
end
xdata   = 0.01:0.01:0.800;
centers = xdata;
fT   = parameterizedData(T,xdata,centers);
t    = fT(xdata);
Udef  = parameterizedData(Td,xdata,centers);
def  = Udef(xdata);
Urb  = parameterizedData(Tr,xdata,centers);
rb   = Urb(xdata);
fK   = parameterizedData(Kcoase,xdata,centers);
Kcoarse = fK(xdata);

% for i=1:length(xdata)
%     error(i) = norm(def(:,i)+rb(:,i)-T(:,i))/norm(T(:,i));
%     error2(i) = norm(t(:,i)-T(:,i))/norm(T(:,i));
% end

%load('data_0.100.mat')
% EIFEoper.U = reshape(def,[],8) + reshape(rb',[],8);
% kcoarse = EIFEoper.U'*EIFEoper.Kfine*EIFEoper.U;
EIFEoper.Kcoarse = fK;
EIFEoper.Udef = Udef;
EIFEoper.Urb  = Urb;

filePath = './EPFL/parametrizedEIFEM.mat';
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

function f = parameterizedData(var,xdata,centers)
deim    = DEIM(var);

coeff   = deim.basis(deim.indices,:)\var(deim.indices,:);
rbf       = RBF(coeff',xdata,centers);
 f = @ (r) deim.basis*rbf.evaluate(r) ;
%f = @ (r) reshape( deim.basis*rbf.evaluate(r),[],8);

end

% If you want to access specific variables in allData, use:
% allData{1}, allData{2}, ..., allData{k}
