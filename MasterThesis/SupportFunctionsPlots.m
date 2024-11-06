
clear;
close all;

n1   = 200;
n2   = 200;
mesh = TriangleMesh(1,1,n1,n2);

x1 = @(x) x(1,:,:)-0.5;
x2 = @(x) x(2,:,:)-0.5;
epsilon = mesh.computeMeanCellSize();
Smagnit = @(x) sqrt(x1(x).^2+x2(x).^2);
% Circle
% h       = @(x) epsilon*sqrt(x1(x).^2+x2(x).^2);
% sig2    = AnalyticalFunction.create(@(x) h(x).^2,1,mesh);
% sig2P1  = sig2.project('P1');
% sig2P1.plot();
% a = gcf().findobj();
% a(4).EdgeColor = 'none';
% sig2P1.print('CircleSupport')

%Ellipse

% theta    = 90;
% alpha    = 0.8;
% A11 = cosd(theta)^2 + sind(theta)^2 / alpha^2;
% A12 = (cosd(theta) * sind(theta) * (alpha^2 - 1)) / alpha^2;
% A21 = (cosd(theta) * sind(theta) * (alpha^2 - 1)) / alpha^2;
% A22 = sind(theta)^2 + cosd(theta)^2 / alpha^2;
% M = [A11, A12; A21, A22];
% A = inv(M);

% h = @(x) epsilon * sqrt(pagemtimes(pagetranspose([x1(x); x2(x)]), pagemtimes(A, [x1(x); x2(x)])));
% h = @(x) epsilon * sqrt(pagemtimes(pagetranspose([(x1(x)./Smagnit(x); (x2(x)./Smagnit(x))]), pagemtimes(A, [(x1(x)./Smagnit(x); (x2(x)./Smagnit(x))])));
% sig2    = AnalyticalFunction.create(@(x) h(x).^2,1,mesh);
% sig2P1D  = sig2.project('P1D');
% sig2P1  = sig2.project('P1');     
% sig2P1D.plotContour();
% a = gcf().findobj();
% a(4).EdgeColor = 'none';
% sig2P1.print('EllipseSupportp2')


%Segment

alpha = 10;
beta = 2;
h = @(x) alpha * (max(0, (x1(x)) * 1 + (x2(x)) * 1)./Smagnit(x)) ...
              - beta * (min(0, (x1(x)) * 1 + (x2(x)) * 1)./Smagnit(x));
s2     = @(x) Smagnit(x)  ./  (h(x));
% sig2    = AnalyticalFunction.create(@(x) h(x),1,mesh);
sig2   = AnalyticalFunction.create(@(x) s2(x),1,mesh);
sig2P1  = sig2.project('P1D');
sig2P1.plotContour();
% a = gcf().findobj();
% a(4).EdgeColor = 'none';
% sig2P1.print('SegmentSupportalph10beta1')


%droplet 1
%alpha = 5.73;%acosd(1/45);
%h = @(x) epsilon* max(sqrt(x1(x).^2+x2(x).^2), (alpha*x1(x) * 0 + alpha*x2(x) * 1));

%sig2    = AnalyticalFunction.create(@(x) h(x),1,mesh);
%sig2P1  = sig2.project('P1D');

 %sig2P1.plotContour();
%a = gcf().findobj();
%a(4).EdgeColor = 'none';
%sig2P1.print('dropletSupport');



% droplet 2
% alpha  = sqrt(2);
% kx     = 0;
% ky     = 1;
% h1     = @(x) max(  1  ,   (alpha*x1(x) * kx + alpha*x2(x) * ky)./sqrt(x1(x).^2+x2(x).^2)  );
% s2     = @(x) sqrt(x1(x).^2+x2(x).^2)  ./  (h1(x).^2);
% sig2   = AnalyticalFunction.create(@(x) s2(x),1,mesh);
% sig2P1 = sig2.project('P1D');
% sig2P1.plotContour();
% a = gcf().findobj();
% a(4).EdgeColor = 'none';