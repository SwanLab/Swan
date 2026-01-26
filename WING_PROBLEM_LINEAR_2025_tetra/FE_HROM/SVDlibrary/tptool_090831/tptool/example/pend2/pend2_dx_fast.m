function dx = pend2_dx(u, x)

m_k = 1; % kg
m_1 = 0.3; % kg
m_2 = 0.2; % kg
g = 9.81; % m/s^2
l_1 = 0.6; % m
l_2 = 0.2; % m

a_1 = x(2);
a_2 = x(3);
da_1= x(5);
da_2= x(6);

E22 = [...
    m_k+m_1+m_2,        m_1*l_1*cos(a_1),   m_2*l_2*cos(a_2);...
    m_1*l_1*cos(a_1),   4/3*m_1*l_1^2,      0;...
    m_2*l_2*cos(a_2),   0,                  4/3*m_2*l_2^2];

E = [...
    eye(3, 3)   zeros(3, 3);...
    zeros(3, 3) E22];

% sin(a)/a
if abs(a_1) < 1e-5
	sa_1 = 1 - a_1^2/6;
else
	sa_1 = sin(a_1)/a_1;
end
if abs(a_2) < 1e-5
	sa_2 = 1 - a_2^2/6;
else
	sa_2 = sin(a_2)/a_2;
end

A21 = [...
    0,  0,                      0;...
    0,  m_1*l_1*g*sa_1, 0;...
    0,  0,                      m_2*l_2*g*sa_2];

% with friction b_k, b_1, b_2
%A22 = [-b_k, m_1*l_1*da_1*sin(a_1), m_2*l_2*da_2*sin(a_2); 0, -b_1, 0; 0, 0, -b_2];
A22 = [...
    0, m_1*l_1*da_1*sin(a_1),   m_2*l_2*da_2*sin(a_2);...
    0, 0,                       0;...
    0, 0,                       0];

A = [ ...
    zeros(3, 3) eye(3, 3); ...
    A21         A22];

B = [0 0 0 1 0 0]';

% E dx = A x + B u
iE = pinv(E);
dx = iE*A*x + iE*B*u;

