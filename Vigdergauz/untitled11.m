theta = sym('theta','real');
phi = sym('phi','real');
ax = sym('ax','real');
ay = sym('ay','real');

eq(1) = theta - (ax/cx + ay/cy - 1);
eq(2) = tan(phi) - (ax/ay);

x = solve(eq,[ax ay]);
pretty(simplify(x.ax))
pretty(simplify(x.ay))


phi = linspace(-pi/2,pi/2,100);
theta = 0.2;
axv = (tan(phi)*(theta + 1))./(tan(phi) + 1);
ayv = (theta + 1)./(tan(phi) + 1);

plot(phi,axv)
hold on
plot(phi,ayv)

