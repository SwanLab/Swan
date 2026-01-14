clear all;

f = @(p)[p(1) p(1)^2]';

% sampling intervals for each parameter
domain = [-2 2];
% grid size: number of grid points for each parameter
gridsize = 19;

% sampling
D = sampling_vec(f, domain, gridsize);

% hosvd
[S U] = hosvd(D, [1 0]);
U = U{1};

[W V] = genhull(U,'box');

Uhat = -ones(size(W));
Uhat(1,1) = 0;
Uhat(8,1) = 0.4;
Uhat(10,1) = 0.6;
Uhat(19,1) = 0;
Uhat(5,2) = 0.6;
Uhat(1,3) = 0;
Uhat(19,3) = 1;
Uhat(1,4) = 1;
Uhat(19,4) = 0;
M = hull_manip(W,Uhat);
Wmod = W*M;
Vmod = M\V;
[i,j] = find(Uhat > -1);
uhat = Uhat(Uhat > -1);
n = length(i);
%for k=1:n
%    [uhat(k) Wmod(i(k),j(k))] %#ok<NOPTS>
%end

close all
set(0, 'DefaultFigureWindowStyle', 'docked')
plothull(U);
plothull(W);
plothull(Wmod);
hold on
plot(i,uhat,'kx')
hold off

figure
hold on
plot(D(:,1), D(:,2))
plot(V*S(:,1), V*S(:,2), 'go')
plot(Vmod*S(:,1), Vmod*S(:,2), 'mo')
hold off
