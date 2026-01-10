function dpp = ppDer(pp)
% dpp = ppDer(pp)
%
% Computes the time-derivative of piece-wise polynomial (PP) struct
%
% INPUTS:
%   pp = a PP struct containing a trajectory of interest
% 
% OUTPUTS:
%   dpp = a new PP struct that is the time-derivative of pp
%
% NOTES:
%   --> a pp struct is typically created by matlab functions such as
%   spline, pchip, or pwch and evaluated using ppval.
%   --> Call this function without arguments to run a test case
%

if nargin == 0
    ppDer_test();
    return;
end

n = pp.order;
nRows = size(pp.coefs,1);
dpp.form = pp.form;
dpp.breaks = pp.breaks;
dpp.coefs = zeros(nRows,n-1);
for i=1:n-1
   dpp.coefs(:,i) = (n-i)*pp.coefs(:,i);
end
dpp.pieces = pp.pieces;
dpp.order = pp.order-1;
dpp.dim = pp.dim;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function ppDer_test()
%

nKnot = 10;

tGrid = linspace(0,10,nKnot);
xGrid = [sin(tGrid);cos(tGrid)];
vGrid = [cos(tGrid);-sin(tGrid)];

pp.x = pwch(tGrid,xGrid,vGrid);
pp.v = ppDer(pp.x);

t = linspace(tGrid(1),tGrid(end),200);
x = ppval(pp.x,t);
v = ppval(pp.v,t);

figure(1); clf;

subplot(2,1,1); hold on;
plot(t,x)
plot(tGrid, xGrid,'ko');
xlabel('time')
ylabel('position')

subplot(2,1,2); hold on;
plot(t,v)
xlabel('time')
ylabel('velocity')
plot(tGrid, vGrid,'ko');

end