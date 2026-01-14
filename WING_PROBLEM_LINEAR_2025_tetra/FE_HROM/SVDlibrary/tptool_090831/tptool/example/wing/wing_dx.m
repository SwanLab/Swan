function dx = wing_dx(u, x)

wing_lpv

v = 20; % 20 m/s (should be the same as in wing_control)
p = [v x(2)];

AB = querylpv(LPV, p);

% dx = A x + B u
dx = AB * [x(:); u];
