function untitled4

theta = sym('theta','positive');
r = sym('r','real');
%r = 1;
h = computeH(theta,r);
h = simplify(h);

cx = 1;
cy = 1;

Tx = computeTx(r,h,cx);
Ty = computeTy(r,h,cx);

ax = computeAx(Tx,Ty,cx);
ay = computeAy(Tx,Ty,cx);

ax = simplify(ax);
ay = simplify(ay);

res = theta - (ax/cx + ay/cy - 1);
res = simplify(res)

res2 = r - (1-ay)/(1-ax/cx);
res2 = simplify(res2)

end

function h = computeH(theta,r)
n = -(1+r)*theta + sqrt(theta^2*(r-1)^2 + 4*r);
d = 2*r*(1+theta);
h = n/d;
end

function Tx = computeTx(r,h,cx)
n = r*h*(1+h);
d = (1+r*h)*cx;
Tx = n/d;
end

function Ty = computeTy(r,h,cx)
n = cx*h*(1+r*h);
d = (1+h);
Ty = n/d;
end

function ax = computeAx(Tx,Ty,cx)
n = cx-Ty;
d = 1-Tx*Ty;
ax = n/d;
end

function ay = computeAy(Tx,Ty,cx)
n = 1-Tx*cx;
d = 1-Tx*Ty;
ay = n/d;
end