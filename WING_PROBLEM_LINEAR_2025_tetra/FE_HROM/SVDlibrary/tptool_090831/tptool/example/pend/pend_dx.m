function dx = pend_dx(u, x)

M = 1;   % kg
m = 0.3; % kg
g = 9.81;% m s^-2
L = 0.6; % m

a = x(2);
da = x(4);

if a == 0
	sinca = 1;
else
	sinca = sin(a)/a;
end

E = [...
	eye(2,2) zeros(2,2);
	zeros(2,2) [M/m+1, L*cos(a); cos(a), 4/3*L]];
A = [...
	zeros(2,2) eye(2,2);
	[0 0; 0 g*sinca] [0 L*da*sin(a); 0 0]];

B = [0 0 1/m 0]';

% E dx = A x + B u
piE = pinv(E);
dx = piE*A*x + piE*B*u;
