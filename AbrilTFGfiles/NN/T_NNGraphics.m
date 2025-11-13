% PERFORMANCE OF THE DOWNSCALING NN AND GRAPHICS

%clc
%clear


%% Loads the data of the NN

filename1='T1.mat';
filePath1 = fullfile('AbrilTFGfiles', 'NN', filename1);
load(filePath1);

% Mesh load
filename2='UL_r0_0000-20x20.mat';
filePath2 = fullfile('AbrilTFGfiles', 'DataVariables', filename2);
load(filePath2, 'mesh');

real=reshape(real.',2,mesh.nnodes,[]);
real=permute(real, [2 1 3]);
predicted=reshape(predicted.',2,mesh.nnodes,[]);
predicted=permute(predicted, [2 1 3]);
difference=reshape(difference.',2,mesh.nnodes,[]);
difference=permute(difference, [2 1 3]);

idx=1; %radius to visualize 

%% Real
r.mesh=mesh;
r.order='P1';
r.fValues=real(:,1,idx);

r.function=LagrangianFunction(r);
r.function.plot();
fig1=gcf;



%% Predicted
p.mesh=mesh;
p.order='P1';
p.fValues=predicted(:,1,idx);

p.function=LagrangianFunction(p);
p.function.plot();
fig2=gcf;


%% Difference

d.mesh=mesh;
d.order='P1';
d.fValues=difference(:,1,idx);

d.function=LagrangianFunction(d);
d.function.plot();
fig3=gcf;


%% Difference ABS

dA.mesh=mesh;
dA.order='P1';
dA.fValues=abs(difference(:,1,idx));
dA.function=LagrangianFunction(dA);
dA.function.plot();
fig4=gcf;


%% Comparison

figure();
t = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');


% Vista d'adalt
nexttile;
copyobj(allchild(get(fig1,'CurrentAxes')), gca);
title('Real');
nexttile;
copyobj(allchild(get(fig2,'CurrentAxes')), gca);
title('Predicted');
nexttile;
copyobj(allchild(get(fig3,'CurrentAxes')), gca);
title('Difference');
nexttile;
copyobj(allchild(get(fig4,'CurrentAxes')), gca);
title('Difference (ABS)');

% Vista d'abaix
nexttile;
copyobj(allchild(get(fig1,'CurrentAxes')), gca);
title('Real');
view(3);
nexttile;
copyobj(allchild(get(fig2,'CurrentAxes')), gca);
title('Predicted');
view(3);
nexttile;
copyobj(allchild(get(fig3,'CurrentAxes')), gca);
title('Difference');
view(3);
nexttile;
copyobj(allchild(get(fig4,'CurrentAxes')), gca);
title('Difference (ABS)');
view(3);

%% VTU FILE TO COMPARE
fileName = ['Real8'];
r.function.print(fileName,'Paraview');