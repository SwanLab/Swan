function u = pend_control(x)

load('pend_data', 'K', 'U', 'domain');

p = [x(2) x(4)];
u = tpcontroller(p, x, K, U, domain);
