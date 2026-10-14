clc
clear
close all

data = ...
    readmatrix('LatticeNNSensitivityData.csv');

X = data(:,1:3);

C = data(:,4:9);

dCda = data(:,10:15);

dCdb = data(:,16:21);

rho = X(:,1);
a = X(:,2);
b = X(:,3);

normC = vecnorm(C,2,2);

normA = vecnorm(dCda,2,2);

normB = vecnorm(dCdb,2,2);

fprintf('\n')
fprintf('============================================\n')
fprintf('      SENSITIVITY DATASET VALIDATION\n')
fprintf('============================================\n')

fprintf('\nNumber of samples = %d\n',size(data,1))

fprintf('\nNorm statistics\n')

fprintf('\nC\n')
fprintf('min  = %.12e\n',min(normC))
fprintf('mean = %.12e\n',mean(normC))
fprintf('max  = %.12e\n',max(normC))

fprintf('\ndC/da\n')
fprintf('min  = %.12e\n',min(normA))
fprintf('mean = %.12e\n',mean(normA))
fprintf('max  = %.12e\n',max(normA))

fprintf('\ndC/db\n')
fprintf('min  = %.12e\n',min(normB))
fprintf('mean = %.12e\n',mean(normB))
fprintf('max  = %.12e\n',max(normB))

fprintf('\n')
fprintf('============================================\n')
fprintf('           RHO ENDPOINT CHECK\n')
fprintf('============================================\n')

isRho0 = abs(rho) < 1e-12;
isRho1 = abs(rho-1) < 1e-12;

fprintf('\nrho = 0\n')
fprintf('max ||dC/da|| = %.12e\n',...
    max(normA(isRho0)))

fprintf('max ||dC/db|| = %.12e\n',...
    max(normB(isRho0)))

fprintf('\nrho = 1\n')
fprintf('max ||dC/da|| = %.12e\n',...
    max(normA(isRho1)))

fprintf('max ||dC/db|| = %.12e\n',...
    max(normB(isRho1)))

fprintf('\n')
fprintf('============================================\n')
fprintf('             LARGEST dC/da\n')
fprintf('============================================\n')

nShow = 10;

[~,orderA] = sort(normA,'descend');

for i = 1:nShow

    id = orderA(i);

    fprintf('\nSample %d\n',id)

    fprintf(...
        'rho = %.6f, a = %.6f, b = %.6f\n',...
        rho(id),a(id),b(id))

    fprintf(...
        '||dC/da|| = %.12e\n',...
        normA(id))

    disp(dCda(id,:))

end

fprintf('\n')
fprintf('============================================\n')
fprintf('             LARGEST dC/db\n')
fprintf('============================================\n')

[~,orderB] = sort(normB,'descend');

for i = 1:nShow

    id = orderB(i);

    fprintf('\nSample %d\n',id)

    fprintf(...
        'rho = %.6f, a = %.6f, b = %.6f\n',...
        rho(id),a(id),b(id))

    fprintf(...
        '||dC/db|| = %.12e\n',...
        normB(id))

    disp(dCdb(id,:))

end

fprintf('\n')
fprintf('============================================\n')
fprintf('       FD CONVERGENCE AT WORST POINTS\n')
fprintf('============================================\n')

hARange = 0.98;
hBRange = 3-0.3;

hRelative = [...
    5e-3,...
    2.5e-3,...
    1.25e-3];

nWorst = 3;

worstIndices = unique([...
    orderA(1:nWorst);
    orderB(1:nWorst)]);

lattice = TutorialLatticeNN(false);

for iPoint = 1:length(worstIndices)

    id = worstIndices(iPoint);

    x = X(id,:);

    fprintf('\n')
    fprintf('--------------------------------------------\n')
    fprintf('Sample %d\n',id)

    fprintf(...
        'rho = %.6f, a = %.6f, b = %.6f\n',...
        x(1),x(2),x(3))

    fprintf('--------------------------------------------\n')

    dA = zeros(length(hRelative),6);
    dB = zeros(length(hRelative),6);

    for iH = 1:length(hRelative)

        hA = ...
            hRelative(iH)*hARange;

        hB = ...
            hRelative(iH)*hBRange;

        xp = x;
        xm = x;

        xp(2) = x(2)+hA;
        xm(2) = x(2)-hA;

        Cp = ...
            lattice.computeHomogenizationAt(...
            xp(1),xp(2),xp(3));

        Cm = ...
            lattice.computeHomogenizationAt(...
            xm(1),xm(2),xm(3));

        Yp = tensorToVector(Cp);
        Ym = tensorToVector(Cm);

        dA(iH,:) = ...
            (Yp-Ym)/(2*hA);

        xp = x;
        xm = x;

        xp(3) = x(3)+hB;
        xm(3) = x(3)-hB;

        Cp = ...
            lattice.computeHomogenizationAt(...
            xp(1),xp(2),xp(3));

        Cm = ...
            lattice.computeHomogenizationAt(...
            xm(1),xm(2),xm(3));

        Yp = tensorToVector(Cp);
        Ym = tensorToVector(Cm);

        dB(iH,:) = ...
            (Yp-Ym)/(2*hB);

    end

    fprintf('\ndC/da\n')

    for iH = 1:length(hRelative)

        fprintf(...
            'h = %.12e   norm = %.12e\n',...
            hRelative(iH)*hARange,...
            norm(dA(iH,:)))

        disp(dA(iH,:))

    end

    convA12 = ...
        norm(dA(2,:)-dA(1,:))/...
        max(norm(dA(2,:)),eps);

    convA23 = ...
        norm(dA(3,:)-dA(2,:))/...
        max(norm(dA(3,:)),eps);

    fprintf(...
        'Convergence h -> h/2   = %.12e\n',...
        convA12)

    fprintf(...
        'Convergence h/2 -> h/4 = %.12e\n',...
        convA23)

    fprintf('\ndC/db\n')

    for iH = 1:length(hRelative)

        fprintf(...
            'h = %.12e   norm = %.12e\n',...
            hRelative(iH)*hBRange,...
            norm(dB(iH,:)))

        disp(dB(iH,:))

    end

    convB12 = ...
        norm(dB(2,:)-dB(1,:))/...
        max(norm(dB(2,:)),eps);

    convB23 = ...
        norm(dB(3,:)-dB(2,:))/...
        max(norm(dB(3,:)),eps);

    fprintf(...
        'Convergence h -> h/2   = %.12e\n',...
        convB12)

    fprintf(...
        'Convergence h/2 -> h/4 = %.12e\n',...
        convB23)

end

function Y = tensorToVector(C)

    Y = [...
        C(1,1,1,1),...
        C(2,2,2,2),...
        C(1,1,2,2),...
        C(1,1,1,2),...
        C(2,2,1,2),...
        C(1,2,1,2)];

end