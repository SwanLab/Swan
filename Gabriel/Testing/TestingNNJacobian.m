X = [0.40 0.30 1.20;
    0.65 0.50 2.00;
    0.25 0.20 0.70];

[Y,J] = T_NN.computeOutputValuesAndGradient(X);

h = 1e-6;

Jfd = zeros(size(J));

for i = 1:size(X,2)

    Xp = X;
    Xm = X;

    Xp(:,i) = Xp(:,i) + h;
    Xm(:,i) = Xm(:,i) - h;

    Yp = T_NN.computeOutputValues(Xp);
    Ym = T_NN.computeOutputValues(Xm);

    dY = (Yp-Ym)/(2*h);

    Jfd(i,:,:) = permute(dY,[3 2 1]);

end

absError = max(abs(J(:)-Jfd(:)));

relError = ...
    norm(J(:)-Jfd(:))/...
    max(norm(Jfd(:)),eps);

fprintf('\nPhysical Jacobian validation\n')
fprintf('Max absolute error = %.12e\n',absError)
fprintf('Relative error     = %.12e\n',relError)

for i = 1:size(X,2)

    Ji = squeeze(J(i,:,:))';
    Jfdi = squeeze(Jfd(i,:,:))';

    errI = max(abs(Ji(:)-Jfdi(:)));

    fprintf('Feature %d: max error = %.12e\n',i,errI)

end