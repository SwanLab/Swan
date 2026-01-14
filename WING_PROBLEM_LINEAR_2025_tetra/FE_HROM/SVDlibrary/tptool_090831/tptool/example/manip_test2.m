clear all;

f = @(p)[p(2)^3+p(1)^2 p(2)+2*p(1)]';

% sampling intervals for each parameter
domain = [-1 1; -1.5 1.5];
% grid size: number of grid points for each parameter
gridsize = [19 17];

% sampling
D = sampling_vec(f, domain, gridsize);

% shift
%c = [-3 0];%mean_vec(D);
%D = shift_vec(D, -c);

% hosvd
[S U] = hosvd(D, [1 1 0], [0.1 0.1 0]);
U = {U{1} U{2}};

[W V] = genhull(U,'box');
Uhat = -1*ones(size(W{1}));
Uhat(1,1) = 0;
Uhat(8,1) = 0.3;
Uhat(10,1) = 0.5;
Uhat(19,1) = 0;
Uhat(5,2) = 0.5;
Uhat(1,3) = 0;
Uhat(19,3) = 1;
Uhat(1,4) = 1;
Uhat(19,4) = 0;
M = hull_manip(W{1},Uhat);
Wmod = {W{1}*M, W{2}};
Vmod = {M\V{1}, V{2}};
[i,j] = find(Uhat > -1);
uhat = Uhat(Uhat > -1);
n = length(i);
%for k=1:n
%    [uhat(k) Wmod{1}(i(k),j(k))] %#ok<NOPTS>
%end

close all
set(0, 'DefaultFigureWindowStyle', 'docked')
plothull(U);
plothull(W);
plothull(Wmod{1});
hold on
plot(i,uhat,'kx')
hold off

figure
dx = D(:,:,1);
dy = D(:,:,2);
hold on
for i=1:gridsize(2)
	plot(dx(:,i), dy(:,i), 'Color', [(i-1)/gridsize(2) 1-(i-1)/gridsize(2) 0])
end

S = tprod(S, Vmod);
dx = S(:,:,1);
dy = S(:,:,2);
for i=1:size(S,2)
	plot(dx(:,i), dy(:,i), 'o', 'Color', [1-(i-1)/size(S,2) 0 i/size(S,2)])
end
hold off
