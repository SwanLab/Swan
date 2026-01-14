clear all;

f = @(p)[p-p^2+3; 2*p];

% sampling intervals for each parameter
domain = [-3 3];
% grid size: number of grid points for each parameter
gridsize = 19;

% sampling
D = sampling_vec(f, domain, gridsize);

% shift
%c = [-3 0];%mean_vec(D);
%D = shift_vec(D, -c);

% hosvd in the first dimension
[S U] = hosvd(D, [1 0]);
U = U{1};

[Wb Vb] = genhull(U,'box');
[Wn Vn] = genhull(U,'cno');
[Wc Vc] = genhull(U,'close');
[Wi Vi] = genhull(U,'irno');

close all
set(0, 'DefaultFigureWindowStyle', 'docked')
plothull(U);
title('hosvd result (coefficients of the basis vectors)')
plothull(Wb);
title('box result (coefficients of the hull vertices)')
plothull(Wn);
title('cno result (coefficients of the hull vertices)')
plothull(Wc);
title('close result (coefficients of the hull vertices)')
plothull(Wi);
title('irno result (coefficients of the hull vertices)')

figure
plot(D(:,1), D(:,2))
hold on
s = zeros(2*size(S,1),2);
s(1:2:end,:) = S;
plot(s(:,1), s(:,2), 'r-')
plot(Vb*S(:,1), Vb*S(:,2), 'go')
plot(Vn*S(:,1), Vn*S(:,2), 'mo')
plot(Vc*S(:,1), Vc*S(:,2), 'ko')
plot(Vi*S(:,1), Vi*S(:,2), 'ro')
hold off
title('data points and hulls')
legend({'data points [fx(p) fy(p)]',
	'singular vectors weighted by singular values',
	'box hull vertices',
	'cno hull vertices',
	'close hull vertices',
	'irno hull vertices'}, 'Location', 'SouthEast')


