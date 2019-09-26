theta = sym('theta','positive');
q = sym('q','real');

ax = sym('ax','real');
ay = sym('ay','real');

eq(1) = ax - (1 + theta + q)/4;
eq(2) = ay - (1 + theta - q)/4;

x = solve(eq,[theta,q]);
x.theta
x.q
