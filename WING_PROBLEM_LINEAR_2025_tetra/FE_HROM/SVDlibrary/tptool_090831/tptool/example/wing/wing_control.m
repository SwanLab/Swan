function u = wing_control(x)

load('wing_data', 'K', 'U', 'domain');

v = 20; % assume 20m/s (= 72km/h) wind speed
p = [v x(2)];

u = tpcontroller(p, x, K, U, domain);
