clc,clear,close all

H = TutorialLatticeNN(false);

rhoFixed = 0.5;
bFixed   = 1;

a = H.a;

nA = length(a);

mat = zeros(2,2,2,2,nA);

for i = 1:nA
    mat(:,:,:,:,i) = ...
        H.computeHomogenizationAt( ...
        rhoFixed,a(i),bFixed);
end

degree = 8;

[f,df,ddf] = ...
    Homogenization1DFitter.computePolynomial( ...
    degree,a,mat);

homogenization.fun   = f;
homogenization.dfun  = df;
homogenization.ddfun = ddf;

save('Case1LatticeRho05B1A', ...
    'mat','a','rhoFixed','bFixed','homogenization')

aPlot = linspace(min(a),max(a),200);

C1111Fit = f{1,1,1,1}(aPlot);
C2222Fit = f{2,2,2,2}(aPlot);
C1122Fit = f{1,1,2,2}(aPlot);
C1112Fit = f{1,1,1,2}(aPlot);
C2212Fit = f{2,2,1,2}(aPlot);
C1212Fit = f{1,2,1,2}(aPlot);

dC1111Fit = df{1,1,1,1}(aPlot);
dC2222Fit = df{2,2,2,2}(aPlot);
dC1122Fit = df{1,1,2,2}(aPlot);
dC1112Fit = df{1,1,1,2}(aPlot);
dC2212Fit = df{2,2,1,2}(aPlot);
dC1212Fit = df{1,2,1,2}(aPlot);

C1111 = squeeze(mat(1,1,1,1,:));
C2222 = squeeze(mat(2,2,2,2,:));
C1122 = squeeze(mat(1,1,2,2,:));
C1112 = squeeze(mat(1,1,1,2,:));
C2212 = squeeze(mat(2,2,1,2,:));
C1212 = squeeze(mat(1,2,1,2,:));

figure

plot(a,C1111,'o')
hold on
plot(aPlot,C1111Fit)
xlabel('a')
ylabel('C_{1111}^h')
title('\rho = 0.5, b = 1')
grid on

figure

plot(a,C2222,'o')
hold on
plot(aPlot,C2222Fit)
xlabel('a')
ylabel('C_{2222}^h')
title('\rho = 0.5, b = 1')
grid on

figure

plot(a,C1122,'o')
hold on
plot(aPlot,C1122Fit)
xlabel('a')
ylabel('C_{1122}^h')
title('\rho = 0.5, b = 1')
grid on

figure

plot(a,C1112,'o')
hold on
plot(aPlot,C1112Fit)
xlabel('a')
ylabel('C_{1112}^h')
title('\rho = 0.5, b = 1')
grid on

figure

plot(a,C2212,'o')
hold on
plot(aPlot,C2212Fit)
xlabel('a')
ylabel('C_{2212}^h')
title('\rho = 0.5, b = 1')
grid on

figure

plot(a,C1212,'o')
hold on
plot(aPlot,C1212Fit)
xlabel('a')
ylabel('C_{1212}^h')
title('\rho = 0.5, b = 1')
grid on

figure

plot(aPlot,dC1111Fit)
xlabel('a')
ylabel('\partial C_{1111}^h / \partial a')
title('\rho = 0.5, b = 1')
grid on

figure

plot(aPlot,dC2222Fit)
xlabel('a')
ylabel('\partial C_{2222}^h / \partial a')
title('\rho = 0.5, b = 1')
grid on

figure

plot(aPlot,dC1122Fit)
xlabel('a')
ylabel('\partial C_{1122}^h / \partial a')
title('\rho = 0.5, b = 1')
grid on

figure

plot(aPlot,dC1112Fit)
xlabel('a')
ylabel('\partial C_{1112}^h / \partial a')
title('\rho = 0.5, b = 1')
grid on

figure

plot(aPlot,dC2212Fit)
xlabel('a')
ylabel('\partial C_{2212}^h / \partial a')
title('\rho = 0.5, b = 1')
grid on

figure

plot(aPlot,dC1212Fit)
xlabel('a')
ylabel('\partial C_{1212}^h / \partial a')
title('\rho = 0.5, b = 1')
grid on

format long g

fprintf('\n============================================\n')
fprintf('CASE: rho = %.2f, b = %.2f, varying a\n', ...
    rhoFixed,bFixed)
fprintf('============================================\n')

fprintf('\n==============================\n')
fprintf('a VALUES\n')
fprintf('==============================\n')

disp(a(:))

fprintf('\n==============================\n')
fprintf('HOMOGENIZED COMPONENTS\n')
fprintf('==============================\n')

results = [ ...
    a(:), ...
    C1111(:), ...
    C2222(:), ...
    C1122(:), ...
    C1112(:), ...
    C2212(:), ...
    C1212(:)];

disp(array2table(results, ...
    'VariableNames',{ ...
    'a', ...
    'C1111', ...
    'C2222', ...
    'C1122', ...
    'C1112', ...
    'C2212', ...
    'C1212'}))

C1111FitData = f{1,1,1,1}(a);
C2222FitData = f{2,2,2,2}(a);
C1122FitData = f{1,1,2,2}(a);
C1112FitData = f{1,1,1,2}(a);
C2212FitData = f{2,2,1,2}(a);
C1212FitData = f{1,2,1,2}(a);

fprintf('\n==============================\n')
fprintf('FITTING AT HOMOGENIZATION POINTS\n')
fprintf('==============================\n')

fitResults = [ ...
    a(:), ...
    C1111(:),C1111FitData(:), ...
    C2222(:),C2222FitData(:), ...
    C1122(:),C1122FitData(:), ...
    C1112(:),C1112FitData(:), ...
    C2212(:),C2212FitData(:), ...
    C1212(:),C1212FitData(:)];

disp(array2table(fitResults, ...
    'VariableNames',{ ...
    'a', ...
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
    a(:), ...
    errorC1111, ...
    errorC2222, ...
    errorC1122, ...
    errorC1112, ...
    errorC2212, ...
    errorC1212];

disp(array2table(errorResults, ...
    'VariableNames',{ ...
    'a', ...
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

dC1111Data = df{1,1,1,1}(a);
dC2222Data = df{2,2,2,2}(a);
dC1122Data = df{1,1,2,2}(a);
dC1112Data = df{1,1,1,2}(a);
dC2212Data = df{2,2,1,2}(a);
dC1212Data = df{1,2,1,2}(a);

fprintf('\n==============================\n')
fprintf('FIRST DERIVATIVES WITH RESPECT TO a\n')
fprintf('==============================\n')

derivativeResults = [ ...
    a(:), ...
    dC1111Data(:), ...
    dC2222Data(:), ...
    dC1122Data(:), ...
    dC1112Data(:), ...
    dC2212Data(:), ...
    dC1212Data(:)];

disp(array2table(derivativeResults, ...
    'VariableNames',{ ...
    'a', ...
    'dC1111da', ...
    'dC2222da', ...
    'dC1122da', ...
    'dC1112da', ...
    'dC2212da', ...
    'dC1212da'}))

fprintf('\n==============================\n')
fprintf('END POINT CHECK\n')
fprintf('==============================\n')

aMin = min(a);
aMax = max(a);

fprintf('\na = %.4f\n',aMin)
fprintf('C1111 = %.12e\n',f{1,1,1,1}(aMin))
fprintf('C2222 = %.12e\n',f{2,2,2,2}(aMin))
fprintf('C1122 = %.12e\n',f{1,1,2,2}(aMin))
fprintf('C1112 = %.12e\n',f{1,1,1,2}(aMin))
fprintf('C2212 = %.12e\n',f{2,2,1,2}(aMin))
fprintf('C1212 = %.12e\n',f{1,2,1,2}(aMin))

fprintf('\na = %.4f\n',aMax)
fprintf('C1111 = %.12e\n',f{1,1,1,1}(aMax))
fprintf('C2222 = %.12e\n',f{2,2,2,2}(aMax))
fprintf('C1122 = %.12e\n',f{1,1,2,2}(aMax))
fprintf('C1112 = %.12e\n',f{1,1,1,2}(aMax))
fprintf('C2212 = %.12e\n',f{2,2,1,2}(aMax))
fprintf('C1212 = %.12e\n',f{1,2,1,2}(aMax))

aMicro = [0 0.2 0.4 0.6 0.8 0.98];

for i = 1:length(aMicro)

    H.plotMicrostructures( ...
        rhoFixed, ...
        aMicro(i), ...
        bFixed);

end