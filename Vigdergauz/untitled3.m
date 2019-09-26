function untitled3

ax = sym('ax','real');
ay = sym('ay','real');
Tx = sym('Tx','real');
Ty = sym('Ty','real');

cx = 1;
cy = 1;

theta = ax/cx + ay/cy - 1;
r = (1-ay/cy)/(1-ax/cx);
r = simplify(r);

h = computeH(theta,r);
h = simplify(h);
pretty(h)

pretty(simplify(theta))

end

function h = computeH(theta,r)
n = -(1+r)*theta + sqrt(theta^2*(r-1)^2 + 4*r);
d = 2*r*(1+theta);
h = n/d;
end