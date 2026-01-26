clear all;

f = @(p)[p-p^2+3 2*p]';

% sampling intervals for each parameter
domain = [-3 3];
% grid size: number of grid points for each parameter
gridsize = 19;

% sampling
D = sampling_vec(f, domain, gridsize);

% hosvd
[S U] = hosvd(D, [1 0]);
U = U{1};

[W V] = genhull(U,'box');
[Wc3 Vc3] = genhull(U);

close all
set(0, 'DefaultFigureWindowStyle', 'docked')
plothull(U);
plothull(W);
plothull(Wc3);

figure
hold on
plot(U(:,1), U(:,2))
plot(V(:,1), V(:,2), 'ko')
plot(Vc3(:,1), Vc3(:,2), 'ro')
plot(Vc3(:,1), Vc3(:,2), 'r-')
hold off


