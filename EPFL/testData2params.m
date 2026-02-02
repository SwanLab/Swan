close all
clear all
% Specify the directory where the .mat files are located
directory = './EPFL/dataLattice2param'; % Update this path as needed

% Get a list of all .mat files in the directory
files = dir(fullfile(directory, 'data_*.mat'));
i=1;
% Loop through each file and load it
for k = 1:length(files)
    % Get the full path to the file
    filePath = fullfile(files(k).folder, files(k).name);

    % Load the file
    load(filePath);

     if EIFEoper.h1 > 0.1 && EIFEoper.h2 > 0.1

          % Optionally, store each loaded variable into a structure or cell array
    % For example, you could store the variables in a cell array:
    T(:,i) = EIFEoper.U;  % This stores each file's contents in the cell array 'allData'
    Td(:,i) = EIFEoper.Udef(:);
    Tr(:,i) = EIFEoper.Urb(:);
    Kcoase(:,i) = EIFEoper.Kcoarse(:);
    Kfine(:,:,i) = full(EIFEoper.Kfine);
    PhiD(:,i) = EIFEoper.PhiD(:);
    PhiR(:,i) = EIFEoper.PhiR(:);
    U(:,i) = EIFEoper.snapshots(:);
     coord(i,:) = [EIFEoper.h1,EIFEoper.h2 ];
     i=i+1;
    % If you want to work with specific variables in each file, you can access them like:
    % data.variableName
     end

   

    disp(['Loaded: ', files(k).name]);  % Display the file being loaded
end

load('parametrizedEIFEMLagrange_2params.mat')
Tdef = EIFEoper.Udef;
Trb = EIFEoper.Urb;
fT   = EIFEoper.U;
fK   = EIFEoper.Kcoarse;

for i=1:size(coord,1)
    aux = fK(coord(i,:)');
    Kcoarse(:,i) = aux(:);
%     error(i) = norm(def(:,i)+rb(:,i)-T(:,i))/norm(T(:,i));
    % error2(i) = norm(t(:,i)-T(:,i))/norm(T(:,i));
    errorK(i) = norm(Kcoarse(:,i)-Kcoase(:,i))/norm(Kcoase(:,i));
%     energyErr(i) = (t(:,i)-T(:,i))'*Kfine(:,:,i)*(t(:,i)-T(:,i));
end
plot(errorK)
% xlim([min(xdata)-0.01,max(xdata)+0.01])
figure
plot(xdata,errorK)
xlim([min(xdata)-0.01,max(xdata)+0.01])
% load('data_0.100.mat')
% EIFEoper.U = reshape(def,[],8) + reshape(rb',[],8);
% kcoarse = EIFEoper.U'*EIFEoper.Kfine*EIFEoper.U;
EIFEoper.Kcoarse = fK;
EIFEoper.Udef = Tdef;
EIFEoper.Urb  = Trb;

