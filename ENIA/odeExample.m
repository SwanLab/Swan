time = [0 5];
fun = @(t,y) 2*t;
y0 = 0;
[t,y] = ode45(fun, time, y0);
plot(t,y)