theta = sym('theta','real');

ax = sym('ax','real');
ay = sym('ay','real');

r = sym('r','real');
cx = 1;
cy = 1;

%h = sym('h','real');

n = -(1+r)*theta + sqrt(theta^2*(r-1)^2 + 4*r);
d = 2*r*(1+theta);
h = n/d;

n = r*h*(1+h);
d = (1+r*h)*cx;
Tx = n/d;

n = cx*h*(1+r*h);
d = (1+h);
Ty = n/d;

axC = simplify((cx - Ty)/(1-Tx*Ty));
ayC = simplify((1 - Tx*cx)/(1-Tx*Ty));


eq(1) = ax - axC;
eq(2) = ay - ayC ;

sol = solve(eq,[theta,r]); 
sol.theta
sol.r

res2 = simplify((r - (axC - axC*ayC)/(ayC - axC*ayC)))
% 
res = theta - (axC/cx + ayC/cy - 1);
res = simplify(res)