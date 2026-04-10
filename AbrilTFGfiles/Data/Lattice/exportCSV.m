close all
clear all
clc;

% Specify case parameters
p.nelem     =  20;
p.Training  = 'EIFEM';         % 'EIFEM'/'Multiscale'
p.Inclusion = 'Material';         % 'Hole'/'Material'/'HoleRaul'
p.Sampling  = 'Isolated';     % 'Isolated'/'Oversampling'
meshName    = p.nelem+"x"+p.nelem;

% Specify the directory where the .mat files are located

% Get a list of all .mat files in the directory
% files = dir(fullfile("AbrilTFGfiles/Data/",p.Training,p.Inclusion,p.Sampling,meshName, 'r0_*.mat'));
files = dir(fullfile("AbrilTFGfiles/Data/Lattice/Multiscale/30x30/", 't1_*.mat'));
% Loop through each file and load it

for k = 1:1:length(files)
    % Get the full path to the file
    filePath = fullfile(files(k).folder, files(k).name);

    % Load the file
    % load(filePath,'EIFEoper',"mesh");
    load(filePath,'Kcoarse','T',"mesh","tCross","tFrame");

    t(k,:)=[tFrame,tCross];

    K_all(:,:,k)=Kcoarse;

    disp(['Loaded: ', files(k).name]);  % Display the file being loaded
end


%% Reshapes the K data and saves it in a csv file

kdata=zeros(size(K_all,2),36+2);
for n=1:size(K_all,3)
    triangSup=triu(K_all(:,:,n));  %gets the triangular superior matrix
    clear row;
    row=[];
    for i=1:8
        for j=i:8
            row(end+1)=triangSup(i,j);
        end
    end
    kdata(n,:)=[t(n,:),row];
end

kFileName = fullfile('AbrilTFGfiles','Data/Lattice/Multiscale/','dataK.csv');
writematrix(kdata,kFileName);

