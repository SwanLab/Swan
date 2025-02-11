clc; clear; close all

s.dim    = '3D';
s.length    = 1;
s.height = 0.1;
test = PerformanceTest(s);
sol = test.compute(0.05);


%% 

% step = 0.05;
% 
% s.dim    = '2D';
% s.length    = 1;
% s.height = 0.1;
% beam     = CantileverBeamMeshCreator(s);
% Nx = fix(s.length/step + 1);
% Ny = fix(s.height/step + 1);
% % Nx = 5;
% % Ny = 2;
% mesh = beam.create(Nx,Ny);
% mesh.plot();

%%
% 
% s.coord = [0 0 0; %1
%            1 0 0; %2
%            0 1 0; %3
%            1 1 0; %4
%            2 1 0; %5
%            2 0 0]; %6
% s.connec = [1 2 3 4;
%             2 4 5 6];
% mesh = Mesh(s);
% mesh.plot

%% Visualize

% Tn = s.connec;
% x  = s.coord(:,1);
% y  = s.coord(:,2);
% figure()
% hold on
% colormap jet;
% plot(x(Tn)',y(Tn)','--','linewidth',0.5);
