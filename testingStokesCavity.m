%% Stokes: Triangle 2D
% file = 'test2d_stokes_triangle_steady';
% a.fileName = file;
% s = StokesDataContainer(a);
% s.material.nu = 1e-3;
% fem = FEM.create(s);
% fem.computeVariables;
% % fem.printPressure(file);
% fem.printVelocity(file);

%% Stokes: Lid-driven cavity
file = 'test2d_stokes_triangle_transient';
a.fileName = file;
s = StokesDataContainer(a);
s.material.nu = 1;
fem = FEM.create(s);
fem.computeVariables;
% fem.printPressure(file);
fem.printVelocity(file,1);

%% Elastic: Triangle
% fileII = 'test2d_triangle';
% b.fileName = fileII;
% d = FemDataContainer(b);
% femII = FEM.create(d);
% femII.solve();
% femII.print(fileII);

%% Elastic: Cantilever
% fileII = 'cantilever3Dfig';
% b.fileName = fileII;
% d = FemDataContainer(b);
% femII = FEM.create(d);
% femII.solve();
% femII.print(fileII);

%% Performance
% clc; clear; close all;
% load('resultsNew.mat')
% load('resultsOld.mat')
% 
% avgN = zeros(numel(resultsNew),1);
% varN = zeros(numel(resultsNew),1);
% avgO = zeros(numel(resultsNew),1);
% varO = zeros(numel(resultsNew),1);
% 
% nelems = [1000, 4000, 18000, 50000, 74000, 119000, 216000, 338000];
% 
% for i = 1:1:numel(resultsNew)
%     resN = resultsNew{i};
%     resO = resultsOld{i};
%     avgN(i) = mean(resN);
%     varN(i) = var(resN);
%     avgO(i) = mean(resO);
%     varO(i) = var(resO);
% end
% 
% figure(1)
% plot(nelems, avgN, '--','linewidth', 1)
% hold on;
% plot(nelems, avgO, '--','linewidth', 1)
% errorbar(nelems, avgN, varN, 'LineStyle','none', 'Color', 'k','linewidth', 0.5);
% errorbar(nelems, avgO, varO, 'LineStyle','none', 'Color', 'k','linewidth', 0.5);
% xlabel('n_{elem}'), ylabel('time (s)'), legend('New','Old');
% title('Performance in a cantilever test')