function u = tora_control(x)

load('tora_data', 'K', 'U', 'domain');

p = [x(3) x(4)];
u = tpcontroller(p, x, K, U, domain);
