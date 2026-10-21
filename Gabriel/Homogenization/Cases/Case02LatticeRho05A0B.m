clc,clear,close all

H = TutorialLatticeNN(false);

rhoFixed = 0.5;
aFixed   = 0;

b = H.b;

nB = length(b);

mat = zeros(2,2,2,2,nB);

for i = 1:nB
    mat(:,:,:,:,i) = ...
        H.computeHomogenizationAt( ...
        rhoFixed,aFixed,b(i));
end

degree = 8;

[f,df,ddf] = ...
    Homogenization1DFitter.computePolynomial( ...
    degree,b,mat);

homogenization.fun   = f;
homogenization.dfun  = df;
homogenization.ddfun = ddf;

folder = fullfile( ...
    'Gabriel', ...
    'HMVademecum', ...
    'Homogenization');

if ~exist(folder,'dir')
    mkdir(folder)
end

fileName = fullfile( ...
    folder, ...
    'Case2LatticeRho05A0B.mat');

save(fileName, ...
    'mat', ...
    'b', ...
    'rhoFixed', ...
    'aFixed', ...
    'homogenization')

bPlot = linspace(min(b),max(b),200);

C1111Fit = f{1,1,1,1}(bPlot);
C2222Fit = f{2,2,2,2}(bPlot);
C1122Fit = f{1,1,2,2}(bPlot);
C1112Fit = f{1,1,1,2}(bPlot);
C2212Fit = f{2,2,1,2}(bPlot);
C1212Fit = f{1,2,1,2}(bPlot);

dC1111Fit = df{1,1,1,1}(bPlot);
dC2222Fit = df{2,2,2,2}(bPlot);
dC1122Fit = df{1,1,2,2}(bPlot);
dC1112Fit = df{1,1,1,2}(bPlot);
dC2212Fit = df{2,2,1,2}(bPlot);
dC1212Fit = df{1,2,1,2}(bPlot);

C1111 = squeeze(mat(1,1,1,1,:));
C2222 = squeeze(mat(2,2,2,2,:));
C1122 = squeeze(mat(1,1,2,2,:));
C1112 = squeeze(mat(1,1,1,2,:));
C2212 = squeeze(mat(2,2,1,2,:));
C1212 = squeeze(mat(1,2,1,2,:));

figure

plot(b,C1111,'o')
hold on
plot(bPlot,C1111Fit)
xlabel('b')
ylabel('C_{1111}^h')
title('\rho = 0.5, a = 0')
grid on

figure

plot(b,C2222,'o')
hold on
plot(bPlot,C2222Fit)
xlabel('b')
ylabel('C_{2222}^h')
title('\rho = 0.5, a = 0')
grid on

figure

plot(b,C1122,'o')
hold on
plot(bPlot,C1122Fit)
xlabel('b')
ylabel('C_{1122}^h')
title('\rho = 0.5, a = 0')
grid on

figure

plot(b,C1112,'o')
hold on
plot(bPlot,C1112Fit)
xlabel('b')
ylabel('C_{1112}^h')
title('\rho = 0.5, a = 0')
grid on

figure

plot(b,C2212,'o')
hold on
plot(bPlot,C2212Fit)
xlabel('b')
ylabel('C_{2212}^h')
title('\rho = 0.5, a = 0')
grid on

figure

plot(b,C1212,'o')
hold on
plot(bPlot,C1212Fit)
xlabel('b')
ylabel('C_{1212}^h')
title('\rho = 0.5, a = 0')
grid on

figure

plot(bPlot,dC1111Fit)
xlabel('b')
ylabel('\partial C_{1111}^h / \partial b')
title('\rho = 0.5, a = 0')
grid on

figure

plot(bPlot,dC2222Fit)
xlabel('b')
ylabel('\partial C_{2222}^h / \partial b')
title('\rho = 0.5, a = 0')
grid on

figure

plot(bPlot,dC1122Fit)
xlabel('b')
ylabel('\partial C_{1122}^h / \partial b')
title('\rho = 0.5, a = 0')
grid on

figure

plot(bPlot,dC1112Fit)
xlabel('b')
ylabel('\partial C_{1112}^h / \partial b')
title('\rho = 0.5, a = 0')
grid on

figure

plot(bPlot,dC2212Fit)
xlabel('b')
ylabel('\partial C_{2212}^h / \partial b')
title('\rho = 0.5, a = 0')
grid on

figure

plot(bPlot,dC1212Fit)
xlabel('b')
ylabel('\partial C_{1212}^h / \partial b')
title('\rho = 0.5, a = 0')
grid on

format long g

fprintf('\n============================================\n')
fprintf('CASE: rho = %.2f, a = %.2f, varying b\n', ...
    rhoFixed,aFixed)
fprintf('============================================\n')

fprintf('\n==============================\n')
fprintf('b VALUES\n')
fprintf('==============================\n')

disp(b(:))

fprintf('\n==============================\n')
fprintf('HOMOGENIZED COMPONENTS\n')
fprintf('==============================\n')

results = [ ...
    b(:), ...
    C1111(:), ...
    C2222(:), ...
    C1122(:), ...
    C1112(:), ...
    C2212(:), ...
    C1212(:)];

disp(array2table(results, ...
    'VariableNames',{ ...
    'b', ...
    'C1111', ...
    'C2222', ...
    'C1122', ...
    'C1112', ...
    'C2212', ...
    'C1212'}))

C1111FitData = f{1,1,1,1}(b);
C2222FitData = f{2,2,2,2}(b);
C1122FitData = f{1,1,2,2}(b);
C1112FitData = f{1,1,1,2}(b);
C2212FitData = f{2,2,1,2}(b);
C1212FitData = f{1,2,1,2}(b);

fprintf('\n==============================\n')
fprintf('FITTING AT HOMOGENIZATION POINTS\n')
fprintf('==============================\n')

fitResults = [ ...
    b(:), ...
    C1111(:),C1111FitData(:), ...
    C2222(:),C2222FitData(:), ...
    C1122(:),C1122FitData(:), ...
    C1112(:),C1112FitData(:), ...
    C2212(:),C2212FitData(:), ...
    C1212(:),C1212FitData(:)];

disp(array2table(fitResults, ...
    'VariableNames',{ ...
    'b', ...
    'C1111','C1111Fit', ...
    'C2222','C2222Fit', ...
    'C1122','C1122Fit', ...
    'C1112','C1112Fit', ...
    'C2212','C2212Fit', ...
    'C1212','C1212Fit'}))

fprintf('\n==============================\n')
fprintf('ABSOLUTE FITTING ERRORS\n')
fprintf('==============================\n')

errorC1111 = C1111FitData(:)-C1111(:);
errorC2222 = C2222FitData(:)-C2222(:);
errorC1122 = C1122FitData(:)-C1122(:);
errorC1112 = C1112FitData(:)-C1112(:);
errorC2212 = C2212FitData(:)-C2212(:);
errorC1212 = C1212FitData(:)-C1212(:);

errorResults = [ ...
    b(:), ...
    errorC1111, ...
    errorC2222, ...
    errorC1122, ...
    errorC1112, ...
    errorC2212, ...
    errorC1212];

disp(array2table(errorResults, ...
    'VariableNames',{ ...
    'b', ...
    'errorC1111', ...
    'errorC2222', ...
    'errorC1122', ...
    'errorC1112', ...
    'errorC2212', ...
    'errorC1212'}))

fprintf('\nMaximum errors:\n')
fprintf('C1111 = %.6e\n',max(abs(errorC1111)))
fprintf('C2222 = %.6e\n',max(abs(errorC2222)))
fprintf('C1122 = %.6e\n',max(abs(errorC1122)))
fprintf('C1112 = %.6e\n',max(abs(errorC1112)))
fprintf('C2212 = %.6e\n',max(abs(errorC2212)))
fprintf('C1212 = %.6e\n',max(abs(errorC1212)))

dC1111Data = df{1,1,1,1}(b);
dC2222Data = df{2,2,2,2}(b);
dC1122Data = df{1,1,2,2}(b);
dC1112Data = df{1,1,1,2}(b);
dC2212Data = df{2,2,1,2}(b);
dC1212Data = df{1,2,1,2}(b);

fprintf('\n==============================\n')
fprintf('FIRST DERIVATIVES WITH RESPECT TO b\n')
fprintf('==============================\n')

derivativeResults = [ ...
    b(:), ...
    dC1111Data(:), ...
    dC2222Data(:), ...
    dC1122Data(:), ...
    dC1112Data(:), ...
    dC2212Data(:), ...
    dC1212Data(:)];

disp(array2table(derivativeResults, ...
    'VariableNames',{ ...
    'b', ...
    'dC1111db', ...
    'dC2222db', ...
    'dC1122db', ...
    'dC1112db', ...
    'dC2212db', ...
    'dC1212db'}))

fprintf('\n==============================\n')
fprintf('END POINT CHECK\n')
fprintf('==============================\n')

bMin = min(b);
bMax = max(b);

fprintf('\nb = %.4f\n',bMin)
fprintf('C1111 = %.12e\n',f{1,1,1,1}(bMin))
fprintf('C2222 = %.12e\n',f{2,2,2,2}(bMin))
fprintf('C1122 = %.12e\n',f{1,1,2,2}(bMin))
fprintf('C1112 = %.12e\n',f{1,1,1,2}(bMin))
fprintf('C2212 = %.12e\n',f{2,2,1,2}(bMin))
fprintf('C1212 = %.12e\n',f{1,2,1,2}(bMin))

fprintf('\nb = %.4f\n',bMax)
fprintf('C1111 = %.12e\n',f{1,1,1,1}(bMax))
fprintf('C2222 = %.12e\n',f{2,2,2,2}(bMax))
fprintf('C1122 = %.12e\n',f{1,1,2,2}(bMax))
fprintf('C1112 = %.12e\n',f{1,1,1,2}(bMax))
fprintf('C2212 = %.12e\n',f{2,2,1,2}(bMax))
fprintf('C1212 = %.12e\n',f{1,2,1,2}(bMax))

bMicro = [0.3 0.6 1 1.5 2 3];

axesMicro = gobjects(length(bMicro),1);
xLimits   = zeros(length(bMicro),2);
yLimits   = zeros(length(bMicro),2);

for i = 1:length(bMicro)

    H.plotMicrostructures( ...
        rhoFixed, ...
        aFixed, ...
        bMicro(i));

    axesMicro(i) = gca;

    hold(axesMicro(i),'on')

    plot(axesMicro(i),0,0,'k+', ...
        'MarkerSize',12, ...
        'LineWidth',1.5)

    text(axesMicro(i),0,0,'  (0,0)', ...
        'VerticalAlignment','bottom')

    axis(axesMicro(i),'equal')

    xLimits(i,:) = xlim(axesMicro(i));
    yLimits(i,:) = ylim(axesMicro(i));

end

xMin = min([xLimits(:,1);0]);
xMax = max([xLimits(:,2);0]);

yMin = min([yLimits(:,1);0]);
yMax = max([yLimits(:,2);0]);

dx = 0.05*(xMax-xMin);
dy = 0.05*(yMax-yMin);

for i = 1:length(axesMicro)

    xlim(axesMicro(i),[xMin-dx,xMax+dx])
    ylim(axesMicro(i),[yMin-dy,yMax+dy])

end