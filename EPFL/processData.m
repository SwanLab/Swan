close all
clear all
% Specify the directory where the .mat files are located
directory = './EPFL/data'; % Update this path as needed

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
xdata   = 1e-6:0.01:0.801;
% N = 80;
% % Interval bounds
% a = 1e-6;
% b = 0.8;
% % Index vector
% i = 0:N;
% % Cosine spacing formula
% xdata = (a + b)/2 + (b - a)/2 * cos(pi * (1 - i / N));
centers = xdata;
[fT,deim]   = parameterizedData(T,xdata,centers);
t    = fT(xdata);
Tdef  = parameterizedData(Td,xdata,centers);
def  = Tdef(xdata);
Trb  = parameterizedData(Tr,xdata,centers);
rb   = Trb(xdata);
fK   = parameterizedData(Kcoase,xdata,centers);
Kcoarse = fK(xdata);

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
EIFEoper.deim    = deim;
filePath = './EPFL/parametrizedEIFEM.mat';
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

load('./EPFL/data2/data_0.79723.mat')
% load('./EPFL/test/data_0.745.mat')
load('./EPFL/data/data_0.100.mat')
EIFEoper.Kcoarse = @(r) EIFEoper.Kcoarse;
EIFEoper.Udef = @(r) EIFEoper.Udef;
EIFEoper.Urb  = @(r) EIFEoper.Urb;

filePath = './EPFL/dataEIFEM.mat';
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

% If you want to access specific variables in allData, use:
% allData{1}, allData{2}, ..., allData{k}
