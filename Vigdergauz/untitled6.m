theta = sym('theta','positive');
r = sym('r','real');

ax = sym('ax','real');
ay = sym('ay','real');

cx = 1;
cy = 1;

eq(1) = theta - (ax/cx + ay/cy - 1);
eq(2) = r - ax/ay*(1-ay)/(1-ax);

x = solve(eq,[ax ay]);
pretty(simplify(x.ax))
pretty(simplify(x.ay))