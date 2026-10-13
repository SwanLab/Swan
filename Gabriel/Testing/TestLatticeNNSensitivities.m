clc
clear
close all

load('LatticeNN.mat','T_NN')

trainDataset = ...
    readmatrix('LatticeNNTrainData.csv');

Xdata = trainDataset(:,1:3);

xMin = min(Xdata,[],1);
xMax = max(Xdata,[],1);
xRange = xMax-xMin;

samplePosition = [...
    0.37 0.43 0.31;
    0.63 0.57 0.74;
    0.48 0.68 0.52];

Xtest = ...
    xMin + samplePosition.*xRange;

hRelative = [2e-2 1e-2 5e-3];

nPoints = size(Xtest,1);
nFeatures = size(Xtest,2);
nOutputs = 6;
nH = length(hRelative);

lattice = TutorialLatticeNN(false);

[Ynn,Jnn] = ...
    T_NN.computeOutputValuesAndGradient(Xtest);

Yfem = zeros(nPoints,nOutputs);

Jfd = zeros(...
    nFeatures,nOutputs,nPoints,nH);

for iPoint = 1:nPoints

    x = Xtest(iPoint,:);

    C0 = lattice.computeHomogenizationAt(...
        x(1),x(2),x(3));

    Yfem(iPoint,:) = ...
        tensorToVector(C0);

    for iFeature = 1:nFeatures

        for iH = 1:nH

            h = ...
                hRelative(iH)*xRange(iFeature);

            xp = x;
            xm = x;

            xp(iFeature) = ...
                xp(iFeature)+h;

            xm(iFeature) = ...
                xm(iFeature)-h;

            Cp = ...
                lattice.computeHomogenizationAt(...
                xp(1),xp(2),xp(3));

            Cm = ...
                lattice.computeHomogenizationAt(...
                xm(1),xm(2),xm(3));

            Yp = tensorToVector(Cp);
            Ym = tensorToVector(Cm);

            Jfd(iFeature,:,iPoint,iH) = ...
                (Yp-Ym)/(2*h);

        end

    end

end

fprintf('\n')
fprintf('============================================\n')
fprintf('        NN VS FEM - TENSOR VALUES\n')
fprintf('============================================\n')

for iPoint = 1:nPoints

    errorC = ...
        norm(Ynn(iPoint,:)-Yfem(iPoint,:))/...
        max(norm(Yfem(iPoint,:)),eps);

    fprintf('\nPoint %d\n',iPoint)

    fprintf(...
        'rho = %.6f, a = %.6f, b = %.6f\n',...
        Xtest(iPoint,1),...
        Xtest(iPoint,2),...
        Xtest(iPoint,3))

    fprintf(...
        'Relative tensor error = %.12e\n',...
        errorC)

end

featureNames = {'rho','a','b'};

fprintf('\n')
fprintf('============================================\n')
fprintf('       NN VS FEM - SENSITIVITIES\n')
fprintf('============================================\n')

for iPoint = 1:nPoints

    fprintf('\n')
    fprintf('--------------------------------------------\n')
    fprintf('Point %d\n',iPoint)

    fprintf(...
        'rho = %.6f, a = %.6f, b = %.6f\n',...
        Xtest(iPoint,1),...
        Xtest(iPoint,2),...
        Xtest(iPoint,3))

    fprintf('--------------------------------------------\n')

    for iFeature = 1:nFeatures

        dNN = reshape(...
            Jnn(iFeature,:,iPoint),...
            1,nOutputs);

        fprintf(...
            '\nDerivative with respect to %s\n',...
            featureNames{iFeature})

        fprintf('\nNN derivative\n')
        disp(dNN)

        for iH = 1:nH

            h = ...
                hRelative(iH)*...
                xRange(iFeature);

            dFD = reshape(...
                Jfd(iFeature,:,iPoint,iH),...
                1,nOutputs);

            absError = ...
                norm(dNN-dFD);

            relError = ...
                absError/...
                max(norm(dFD),eps);

            fprintf(...
                'h = %.6e   Relative error = %.6e\n',...
                h,relError)

        end

        dFD1 = reshape(...
            Jfd(iFeature,:,iPoint,nH-1),...
            1,nOutputs);

        dFD2 = reshape(...
            Jfd(iFeature,:,iPoint,nH),...
            1,nOutputs);

        fdConvergence = ...
            norm(dFD2-dFD1)/...
            max(norm(dFD2),eps);

        fprintf(...
            'FD convergence last two h = %.6e\n',...
            fdConvergence)

    end

end

fprintf('\n')
fprintf('============================================\n')
fprintf('        FINAL DERIVATIVE COMPARISON\n')
fprintf('============================================\n')

for iPoint = 1:nPoints

    fprintf('\n')
    fprintf('--------------------------------------------\n')
    fprintf('Point %d\n',iPoint)

    fprintf(...
        'rho = %.6f, a = %.6f, b = %.6f\n',...
        Xtest(iPoint,1),...
        Xtest(iPoint,2),...
        Xtest(iPoint,3))

    fprintf('--------------------------------------------\n')

    for iFeature = 1:nFeatures

        dNN = reshape(...
            Jnn(iFeature,:,iPoint),...
            1,nOutputs);

        dFEM = reshape(...
            Jfd(iFeature,:,iPoint,end),...
            1,nOutputs);

        absError = ...
            norm(dNN-dFEM);

        relError = ...
            absError/...
            max(norm(dFEM),eps);

        fprintf(...
            '\n%s\n',...
            featureNames{iFeature})

        fprintf(...
            '||dC_NN||  = %.12e\n',...
            norm(dNN))

        fprintf(...
            '||dC_FEM|| = %.12e\n',...
            norm(dFEM))

        fprintf(...
            'Abs error  = %.12e\n',...
            absError)

        fprintf(...
            'Rel error  = %.12e\n',...
            relError)

        disp('NN')
        disp(dNN)

        disp('FEM')
        disp(dFEM)

    end

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