
% Fast computation of the vector q defined by 
% q_k = \sum_{l = 1}^{N_y} G(X_k - Y_l) f_l, k = 1 .. Nx
% Using the "Efficient Bessel Decomposition"
% Code developped by Martin Averseng
% See also SCSD method by François Alouges and Matthieu Aussal. 


addpath(genpath(pwd)); % Add folders of the toolbox to the path. 
clear all;%#ok
close all;
clc;


Nx = 10^4;
Ny = 10^4;
fprintf('10^%s x 10^%s cloud of points \n \n',num2str(round(log(Nx)/log(10))),num2str(round(log(Ny)/log(10))))
% Data points
X = uniformCircle([0,0],1,Nx);
Y = X; %uniformCircle([0,0],1,Ny)+0.5; %uniformDisk([0.2,0],1,Ny);
f = rand(size(Y,1),1); % Vector f

% Parameter for rescaling
rMax = rMaxCalc(X,Y);


% Kernel choice:

% G = LogKernel; % G(x) = log(x)

G = Y0Kernel(0.3); % G(x) = Y0(0.3*x) => Bessel decomposition with Robin 
% conditions. 

% G = 3*Y0Kernel(4); % G(x) = Y0(2*x) => Method of rescaling to a root of Y0

% G = H0Kernel(1000); % G(x) = Y0(1000*x) => Selects frequencies near 0 and 1000


% G = Kernel(@(r)(exp(-r.^2)),@(r)(-2*r.*exp(-r.^2))); % Arbitrary (smooth)
% kernel

% G = Kernel(@(r)(1./r.^2 ),@(r)(-2./r.^3)); % Arbitrary (singular) kernel


% Choice of the cutoff parameter. 
lambda = 3;
a = lambda/sqrt(sqrt(Nx*Ny)); %this value is of the order of the optimal 
% value for data uniformly distributed in a disk. Choose lambda by trial
% and error to minimize the online time.

tol = 1e-3; % input tolerance

% Offline computations.
gradOpt = true;
timeOff = tic;
[dxGxy,dyGxy,rq,loc] = offline_dEBD(G,X,Y,a,tol); 
timeOff = toc(timeOff);
fprintf('Offline computations performed in \n %s seconds \n',num2str(timeOff))
% show the radial quadrature : 
% rq.showDer;

% Online procedure.
timeOn = tic;
qx = dxGxy(f);
qy = dyGxy(f);
timeOn = toc(timeOn);
fprintf('Online products computed in \n %s seconds \n',num2str(timeOn))

% Error on first entry of q 
dist = sqrt((X(1,1) - Y(:,1)).^2 + (X(1,2) - Y(:,2)).^2);
distInv = 1./dist;
distInv(dist<1e-12) = 0;
Y_X_x = (Y(:,1) - X(1,1));
Y_X_y = (Y(:,2) - X(1,2));
qvalx = sum((G.evalDer(dist).*Y_X_x.*distInv).*f);
qvaly = sum((G.evalDer(dist).*Y_X_y.*distInv).*f);
disp('error on first entry');
disp(abs(qvalx - qx(1))/(norm(qx,1)));
disp(abs(qvaly - qy(1))/(norm(qy,1)));

% Error when f = [1 0 0 ... 0]
dist = sqrt((Y(1,1) - X(:,1)).^2 + (Y(1,2) - X(:,2)).^2);
distInv = 1./dist;
distInv(dist < 1e-12) = 0;
Y_X_x = (Y(1,1) - X(:,1));
Y_X_y = (Y(1,2) - X(:,2));
qvalx = (G.evalDer(dist).*Y_X_x.*distInv);
qvaly = (G.evalDer(dist).*Y_X_y.*distInv);
f = [1; zeros(size(Y,1)-1,1)];

qx = dxGxy(f);
qy = dyGxy(f); 

fprintf('Linf error for f = [1 0 0 ... 0] \n (effective error / target accuracy) \n');
fprintf('%s / %s \n\n',num2str(max(abs(qvalx - qx))),num2str(tol))
fprintf('%s / %s \n\n',num2str(max(abs(qvaly - qy))),num2str(tol))

disp('Done');

