

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale' 'EIFisol'
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'
p.nelem      = 30;
meshName     = p.nelem+"x"+p.nelem;

%NN
filePath = fullfile("AbrilTFGfiles","Data/Lattice/",p.Training,"T_NN.mat");
load(filePath);
pol_deg1=pol_deg;

% Dataset
directory= fullfile("AbrilTFGfiles/Data/Lattice/",p.Training,meshName);
files = dir(fullfile(directory, 't1_*.mat'));
i=1;
for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    data=load(filePath);
    training.T(:,i) = data.T(:);  % This stores each file's contents in the cell array 'allData'
    training.mesh{i}=data.mesh;
    training.vec(i,:)=[data.tFrame,data.tCross];
    i=i+1;
end


mesh=training.mesh{1};
training.T1= zeros(mesh.nnodes*mesh.ndim*8,length(files));

for i=1:length(files)
    aux=computeT_NN(training.mesh{i},training.vec(i,:),T_NN,pol_deg1);
    training.T1(:,i)=aux(:);
end

err=vecnorm(abs(training.T-training.T1))/norm(training.T);

figure
scatter3(training.vec(:,1),training.vec(:,2),err);
title("Training Error vs r");
xlabel('t1');
ylabel('t2');
ylim([0 0.035])

tri=delaunay(training.vec(:,1),training.vec(:,2));
trisurf(tri,training.vec(:,1),training.vec(:,2),err);
