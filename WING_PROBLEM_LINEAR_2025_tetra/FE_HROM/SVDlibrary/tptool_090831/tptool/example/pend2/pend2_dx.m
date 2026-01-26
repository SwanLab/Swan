function dx = pend2_dx(u, x)

pend2_lpv

p = [x(2) x(3) x(5) x(6)];
AB = querylpv(LPV, p);

% dx = A x + B u
dx = AB*[x; u];
