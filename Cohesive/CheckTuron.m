Df = 1;
Dc = 0.1;
lambda = linspace(0,2,100);
L = @(l) Df * (l    - Dc)./ (l)/ (Df-Dc);

% plot(lambda,L(lambda));

%lambda = 0.8
%d = 0.9722
