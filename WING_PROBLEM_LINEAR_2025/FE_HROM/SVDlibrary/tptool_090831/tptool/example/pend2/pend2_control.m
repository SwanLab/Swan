function u = pend2_control(x)

load('pend2_data', 'K', 'U', 'domain');

p = [x(2) x(3) x(5) x(6)];
u = tpcontroller(p, x, K, U, domain);
