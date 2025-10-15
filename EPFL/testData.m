close all
clear all
% Specify the directory where the .mat files are located
directory = '/home/raul/Documents/GitHub/EPFL/test'; % Update this path as needed

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
% T=T./vecnorm(T);
load('parametrizedEIFEM.mat')
Tdef = EIFEoper.Udef;
Trb = EIFEoper.Urb;
fT   = EIFEoper.U;

coeff  = EIFEoper.deim.basis(EIFEoper.deim.indices,:)\T(EIFEoper.deim.indices,:);
t= EIFEoper.deim.basis*coeff;
xdata   = 0.005:0.01:0.75;
% t    = fT(xdata);
def  = Tdef(xdata);
rb   = Trb(xdata);
Kcoarse = EIFEoper.Kcoarse(xdata);
% fK   = parameterizedData(Kcoase,xdata,centers);
% Kcoarse = fK(xdata);

for i=1:length(xdata)
    error(i) = norm(def(:,i)+rb(:,i)-T(:,i))/norm(T(:,i));
    error2(i) = norm(t(:,i)-T(:,i))/norm(T(:,i));
    errorK(i) = norm(Kcoarse(:,i)-Kcoase(:,i))/norm(Kcoase(:,i));
end

% load('data_0.100.mat')
% EIFEoper.U = reshape(def,[],8) + reshape(rb',[],8);
% kcoarse = EIFEoper.U'*EIFEoper.Kfine*EIFEoper.U;
EIFEoper.Kcoarse = fK;
EIFEoper.Udef = Tdef;
EIFEoper.Urb  = Trb;

filePath = '/home/raul/Documents/GitHub/EPFL/parametrizedEIFEM.mat';
save(filePath,'EIFEoper')

% load('data_0.800.mat')
% def  = Tdef(0.1);
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

load('data_0.100.mat')
EIFEoper.Kcoarse = @(r) EIFEoper.Kcoarse;
EIFEoper.Udef = @(r) EIFEoper.Udef;
EIFEoper.Urb  = @(r) EIFEoper.Urb;

filePath = '/home/raul/Documents/GitHub/EPFL/dataEIFEM.mat';
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
% f = @ (r) deim.basis*rbf.evaluate(r) ;
f = @ (r) reshape( deim.basis*rbf.evaluate(r),[],8);

end

% If you want to access specific variables in allData, use:
% allData{1}, allData{2}, ..., allData{k}
