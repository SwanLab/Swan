k = sym('k','real');
m = sym('m','real');
N = sym('N','real');
E  = sym('E','real');
nu = sym('nu','real');

eq(1) = k - E./(N*(1-(N-1)*nu));
eq(2) = m - E./(2*(1+nu));

s = solve(eq,[E,nu]);
pretty(simplify(s.E))
pretty(simplify(s.nu))
