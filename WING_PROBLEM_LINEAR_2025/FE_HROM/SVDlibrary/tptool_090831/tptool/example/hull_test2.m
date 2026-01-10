clear all;

f = @(p)[p^3 p]';

% sampling intervals for each parameter
domain = [-2 2];
% grid size: number of grid points for each parameter
gridsize = 17;

% sampling
D = sampling_vec(f, domain, gridsize);

% shift
%...

% hosvd
[S U] = hosvd(D, [1 0]);
U = U{1};

Umax = max(U);
Umin = min(U);

[Wb Vb] = genhull(U,'box');
[Wc Vc] = genhull(U,'cno');

close all
set(0, 'DefaultFigureWindowStyle', 'docked')
plothull(U);
title('hosvd result (coefficients of the basis vectors)')
plothull(Wb);
title('box result (coefficients of the hull vertices)')
plothull(Wc);
title('cno result (coefficients of the hull vertices)')

figure
plot(D(:,1), D(:,2))
hold on
s = zeros(2*size(S,1),2);
s(1:2:end,:) = S;
%plot(s(:,1), s(:,2), 'r-')
plot(Vb*S(:,1), Vb*S(:,2), 'go')
plot(Vc*S(:,1), Vc*S(:,2), 'mo')
hold off
title('data points and hulls')
legend({'data points [fx(p) fy(p)]', 'singular vectors weighted by singular values', 'box hull vertices', 'cno hull vertices'}, 'Location', 'SouthEast')

