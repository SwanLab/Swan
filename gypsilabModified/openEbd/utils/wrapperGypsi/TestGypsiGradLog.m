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

kernel_string = 'grady[log(r)]';
tOff = tic;
[dGxy,loc,K] = MVproduct(kernel_string,X,Y,tol,'lambda',5);
tOff = toc(tOff);
fprintf('Offline computations : %s s\n\n',num2str(tOff));

% Online procedure.
timeOn = tic;
qx = dGxy{1}(f);
qy = dGxy{2}(f);
timeOn = toc(timeOn);
fprintf('Online computations : %s seconds \n',num2str(timeOn))

% Error on first entry of q 
dist = sqrt((X(1,1) - Y(:,1)).^2 + (X(1,2) - Y(:,2)).^2);
distInv = 1./dist;
distInv(dist<1e-12) = 0;
Y_X_x = (Y(:,1) - X(1,1));
Y_X_y = (Y(:,2) - X(1,2));
qvalx = sum((K.evalDer(dist).*Y_X_x.*distInv).*f);
qvaly = sum((K.evalDer(dist).*Y_X_y.*distInv).*f);
disp('error on first entry');
disp(abs(qvalx - qx(1))/(norm(qx,1)));
disp(abs(qvaly - qy(1))/(norm(qy,1)));

% Error when f = [1 0 0 ... 0]
dist = sqrt((Y(1,1) - X(:,1)).^2 + (Y(1,2) - X(:,2)).^2);
distInv = 1./dist;
distInv(dist < 1e-12) = 0;
Y_X_x = (Y(1,1) - X(:,1));
Y_X_y = (Y(1,2) - X(:,2));
qvalx = (K.evalDer(dist).*Y_X_x.*distInv);
qvaly = (K.evalDer(dist).*Y_X_y.*distInv);
f = [1; zeros(size(Y,1)-1,1)];
qx = dGxy{1}(f);
qy = dGxy{2}(f); 
fprintf('Linf error for f = [1 0 0 ... 0] \n (effective error / target accuracy) \n');
fprintf('%s / %s \n\n',num2str(max(abs(qvalx - qx))),num2str(tol))
fprintf('%s / %s \n\n',num2str(max(abs(qvaly - qy))),num2str(tol))
disp('Michto Gypsilab! ;)');

