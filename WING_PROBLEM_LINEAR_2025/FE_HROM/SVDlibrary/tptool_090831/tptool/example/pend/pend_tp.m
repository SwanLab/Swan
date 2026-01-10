% parameters: angle, dangle/dt
domain = [-pi/12 pi/12; -pi pi];
gridsize = [45; 15];

M = 1.0; % kg
m = 0.3; % kg
g = 9.8; % m s^-2
L = 0.6; % m

% LPV model
F = @(x)x(2)*sin(x(1));
G = @(x)g*sinc(x(1)/pi);
H = @(x)-4/3*(M/m+1)+cos(x(1))^2;
lpv = {...
	@(x)0, @(x)0,                      @(x)1, @(x)0,                   @(x)0;
	@(x)0, @(x)0,                      @(x)0, @(x)1,                   @(x)0;
	@(x)0, @(x)G(x)*cos(x(1))/H(x),    @(x)0, @(x)-4/3*L*F(x)/H(x),    @(x)-4/(3*m*H(x));
	@(x)0, @(x)-(M+m)*G(x)/(m*L*H(x)), @(x)0, @(x)F(x)*cos(x(1))/H(x), @(x)cos(x(1))/(m*L*H(x))};

% parameter dependency
dep = zeros([size(lpv) 2]);
dep(3,2,:) = [1 0];
dep(3,4,:) = [1 1];
dep(3,5,:) = [1 0];
dep(4,2,:) = [1 0];
dep(4,4,:) = [1 1];
dep(4,5,:) = [1 0];

% size of the A matrix
n = 4;

[S U] = tptrans(lpv, dep, domain, gridsize, 'close');

plothull(U);
[maxerr, meanerr] = tperror(lpv, S, U, domain, 100);
disp('max and mean error:'); disp(maxerr); disp(meanerr);

save('pend_data', 'S', 'U', 'n', 'domain', 'gridsize');
