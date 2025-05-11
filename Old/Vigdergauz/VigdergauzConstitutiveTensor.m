function VigdergauzConstitutiveTensor
s1 = sym('s1','real');
s2 = sym('s2','real');
s3 = sym('s3','real');

e1 = sym('e1','real');
e2 = sym('e2','real');
e3 = sym('e3','real');

k1 = sym('k1','real');
k2 = sym('k2','real');

mu1 = sym('mu1','real');
mu2 = sym('mu2','real');

t1 = sym('t1','real');
t2 = sym('t2','real');

A = sym('A','real');

lhs = s1 + s2;
rhs = A*(e1 + e2);
eq(1) = lhs - rhs;

lhs = s2 - s1;
rhs = 2*mu2*(e2 - e1);
eq(2) = lhs - rhs;

lhs = s3;
rhs = 2*mu2*(e3);
eq(3) = lhs - rhs;


x = solve(eq,[s1,s2,s3]);
s(1) = x.s1;
s(2) = x.s2;
s(3) = x.s3;

e(1) = e1;
e(2) = e2;
e(3) = e3;


for i = 1:3
    for j = 1:3
        C(i,j) = diff(s(i),e(j));
    end
end


n = k1*k2+mu2*(k1*t1+k2*t2);
d = mu2+t1*k2+t2*k1;

Ar = n/d;

C = subs(C,A,Ar);
C = subs(C,t1,1-t2);
C = simplify(C)
pretty(C);

pretty(inv(0.5*[t1 + s1, t1-s1;t1-s1,t1+s1]))


end