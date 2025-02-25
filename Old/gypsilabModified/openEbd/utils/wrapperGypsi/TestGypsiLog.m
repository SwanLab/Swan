clc;
clear all;
close all;

disp('hello');
N = 10^4;
tol = 1e-2;
fprintf('10^%s x 10^%s cloud of points \n \n',num2str(round(log(N)/log(10))),num2str(round(log(N)/log(10))))
% Data points
X = uniformCircle([0,0],1,N);
Y = X; %uniformDisk([0.2,0],1,Ny);

f = rand(size(Y,1),1); % Vector f

kernel_string = '[log(r)]';
tOff = tic;
[Gxy,loc,K] = MVproduct(kernel_string,X,Y,tol,'lambda',5);
tOff = toc(tOff);
fprintf('Offline computations : %s s\n\n',num2str(tOff));

% Online procedure.
timeOn = tic;
q = Gxy(f);
timeOn = toc(timeOn);
fprintf('Online computations : %s seconds \n',num2str(timeOn))

% Error on first entry of q 
dist = sqrt((X(1,1) - Y(:,1)).^2 + (X(1,2) - Y(:,2)).^2);
qval = sum(K.eval(dist).*f);
disp('error on first entry');
disp(abs(qval - q(1))/(norm(q,1)));

% Error when f = [1 0 0 ... 0]
dist = sqrt((X(:,1) - Y(1,1)).^2 + (X(:,2) - Y(1,2)).^2);
f = [1; zeros(size(Y,1)-1,1)];
q = Gxy(f);
qval = K.eval(dist);
fprintf('Linf error for f = [1 0 0 ... 0] \n (effective error / target accuracy) \n');
fprintf('%s / %s \n\n',num2str(max(abs(qval - q))),num2str(tol))
disp('Michto Gypsilab! ;)');