function dx = tora_dx(u, x)

tora_lpv

p = [x(3) x(4)];

AB = querylpv(lpv, p);

% dx = A x + B u
dx = AB * [x(:); u];
