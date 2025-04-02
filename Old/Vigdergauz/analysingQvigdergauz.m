Fr1 = sym('Fr1','real');
F1r = sym('F1r','real');
Fr2 = sym('Fr2','real');
F2r = sym('F2r','real');

f1 = sym('f1','real');
f2 = sym('f2','real');


t1 = F1r/Fr1;
t2 = F2r/Fr2;

s1 = f1/Fr1;
s2 = f2/Fr2;


% t1 = sym('t1','real');
% t2 = sym('t2','real');
% s1 = sym('s1','real');
% s2 = sym('s2','real');

m1 = (1 - t2)/(1-t1*t2)*s1;
m2 = (1 - t1)/(1-t1*t2)*s2;

volVoid = (t1*(1-t2) + t2*(1-t1))/(1-t1*t2);

cq = volVoid/(m1*m2);
cq = simplify(cq);
pretty(cq)

