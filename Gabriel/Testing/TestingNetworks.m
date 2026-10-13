clc,clear,close all

rng(1)

data.nFeatures = 3;
data.Xtrain = zeros(10,3);
data.nLabels = 6;

s.data = data;
s.hiddenLayers = [8 10 7];
s.HUtype = 'ReLU';
s.OUtype = 'linear';

net = Network(s);

X = [0.2 0.3 1.0;
    0.7 0.5 1.5;
    0.4 0.1 2.0;
    0.8 0.6 0.7];

[Y,J] = net.computeYOutAndJacobian(X);

h = 1e-6;

Jfd = zeros(size(J));

for i = 1:size(X,2)

    Xp = X;
    Xm = X;

    Xp(:,i) = Xp(:,i) + h;
    Xm(:,i) = Xm(:,i) - h;

    Yp = net.computeYOut(Xp);
    Ym = net.computeYOut(Xm);

    dY = (Yp-Ym)/(2*h);

    Jfd(i,:,:) = permute(dY,[3 2 1]);

end

absError = max(abs(J(:)-Jfd(:)));

relError = norm(J(:)-Jfd(:))/max(norm(Jfd(:)),eps);

disp(' ')
disp('Size X')
disp(size(X))

disp('Size Y')
disp(size(Y))

disp('Size Jacobian')
disp(size(J))

fprintf('Max absolute error = %.12e\n',absError)
fprintf('Relative error     = %.12e\n',relError)

for i = 1:size(X,2)

    Ji = squeeze(J(i,:,:))';
    Jfdi = squeeze(Jfd(i,:,:))';

    errI = max(abs(Ji(:)-Jfdi(:)));

    fprintf('Feature %d: max error = %.12e\n',i,errI)

end