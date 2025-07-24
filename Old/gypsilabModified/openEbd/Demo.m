
% Fast computation of the vector q defined by 
% q_k = \sum_{l = 1}^{N_y} G(X_k - Y_l) f_l, k = 1 .. Nx
% Using the "Efficient Bessel Decomposition"
% Code developped by Martin Averseng
% See also SCSD method by François Alouges and Matthieu Aussal. 


addpath(genpath(pwd)); % Add folders of the toolbox to the path. 
clear all;%#ok
close all;
clc;


Nx = 10^3;
Ny = 10^3;
fprintf('10^%s x 10^%s cloud of points \n \n',num2str(round(log(Nx)/log(10))),num2str(round(log(Ny)/log(10))))
% Data points
X = uniformDisk([0,0],1,Nx);
Y = X; %uniformDisk([0.2,0],1,Ny);
f = rand(size(Y,1),1); % Vector f

% Parameter for rescaling
rMax = rMaxCalc(X,Y);


% Kernel choice:

% G = 3*LogKernel; % G(x) = log(x)

% G = Y0Kernel(0.1); % G(x) = Y0(0.1*x) => Bessel decomposition with Robin 
%conditions. 

% G = 3*Y0Kernel(4); % G(x) = Y0(2*x) => Method of rescaling to a root of Y0

% G = H0Kernel(1000); % G(x) = Y0(1000*x) => Selects frequencies near 0 and 1000



G = ThinPlate(10,25); % G(x) = 10*x^2*log(25*x)

% G = Kernel(@(r)(exp(-r.^2)),@(r)(-2*r.*exp(-r.^2))); % Arbitrary (smooth)
% kernel

% G = Kernel(@(r)(1./r.^2 ),@(r)(-2./r.^3)); % Arbitrary (singular) kernel


% Choice of the cutoff parameter. 
lambda = 5;
a = lambda/sqrt(sqrt(Nx*Ny)); %this value is of the order of the optimal 
% value for data uniformly distributed in a disk. Choose lambda by trial
% and error to minimize the online time.

tol = 1e-3; % input tolerance

% Offline computations.
timeOff = tic;
[onlineEBD,rq,loc] = offlineEBD(G,X,Y,a,tol); 
timeOff = toc(timeOff);
fprintf('Offline computations performed in \n %s seconds \n',num2str(timeOff))
% show the radial quadrature : 
rq.show;

% Online procedure.
timeOn = tic;
q = onlineEBD(f);
timeOn = toc(timeOn);
fprintf('Online product computed in \n %s seconds \n',num2str(timeOn))

% Error on first entry of q 
dist = sqrt((X(1,1) - Y(:,1)).^2 + (X(1,2) - Y(:,2)).^2);
qval = sum(G.eval(dist).*f);
disp('error on first entry');
disp(abs(qval - q(1))/(norm(q,1)));

% Error when f = [1 0 0 ... 0]
dist = sqrt((X(:,1) - Y(1,1)).^2 + (X(:,2) - Y(1,2)).^2);
f = [1; zeros(size(Y,1)-1,1)];
q = onlineEBD(f);
qval = G.eval(dist);
fprintf('Linf error for f = [1 0 0 ... 0] \n (effective error / target accuracy) \n');
fprintf('%s / %s \n\n',num2str(max(abs(qval - q))),num2str(tol))
disp('Done');

