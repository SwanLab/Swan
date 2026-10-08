clc,clear,close all

H = TutorialHomogenization();

rho = H.rho;
mat = H.Chomog;

degree = 8;

[f,df,ddf] = DensityHomogenizationFitter.computePolynomial(degree,rho,mat);

homogenization.fun   = f;
homogenization.dfun  = df;
homogenization.ddfun = ddf;

save('SquareDensity','mat','rho','homogenization')

rhoPlot = linspace(0,1,200);

C11Fit = f{1,1,1,1}(rhoPlot);
C22Fit = f{2,2,2,2}(rhoPlot);
C12Fit = f{1,1,2,2}(rhoPlot);
C33Fit = f{1,2,1,2}(rhoPlot);

dC11Fit = df{1,1,1,1}(rhoPlot);
dC22Fit = df{2,2,2,2}(rhoPlot);
dC12Fit = df{1,1,2,2}(rhoPlot);
dC33Fit = df{1,2,1,2}(rhoPlot);

C11 = squeeze(mat(1,1,1,1,:));
C22 = squeeze(mat(2,2,2,2,:));
C12 = squeeze(mat(1,1,2,2,:));
C33 = squeeze(mat(1,2,1,2,:));

figure

plot(rho,C11,'o')
hold on
plot(rhoPlot,C11Fit)
xlabel('\rho')
ylabel('C_{1111}^h')
grid on

figure

plot(rho,C22,'o')
hold on
plot(rhoPlot,C22Fit)
xlabel('\rho')
ylabel('C_{2222}^h')
grid on

figure

plot(rho,C12,'o')
hold on
plot(rhoPlot,C12Fit)
xlabel('\rho')
ylabel('C_{1122}^h')
grid on

figure

plot(rho,C33,'o')
hold on
plot(rhoPlot,C33Fit)
xlabel('\rho')
ylabel('C_{1212}^h')
grid on

figure

plot(rhoPlot,dC11Fit)
xlabel('\rho')
ylabel('\partial C_{1111}^h / \partial \rho')
grid on

figure

plot(rhoPlot,dC22Fit)
xlabel('\rho')
ylabel('\partial C_{2222}^h / \partial \rho')
grid on

figure

plot(rhoPlot,dC12Fit)
xlabel('\rho')
ylabel('\partial C_{1122}^h / \partial \rho')
grid on

figure

plot(rhoPlot,dC33Fit)
xlabel('\rho')
ylabel('\partial C_{1212}^h / \partial \rho')
grid on

format long g

fprintf('\n==============================\n')
fprintf('DENSITY VALUES\n')
fprintf('==============================\n')
disp(rho(:))

fprintf('\n==============================\n')
fprintf('HOMOGENIZED COMPONENTS\n')
fprintf('==============================\n')

results = [rho(:),C11(:),C22(:),C12(:),C33(:)];

disp(array2table(results,...
    'VariableNames',{'rho','C11','C22','C12','C33'}))

C11FitData = f{1,1,1,1}(rho);
C22FitData = f{2,2,2,2}(rho);
C12FitData = f{1,1,2,2}(rho);
C33FitData = f{1,2,1,2}(rho);

fprintf('\n==============================\n')
fprintf('FITTING AT HOMOGENIZATION POINTS\n')
fprintf('==============================\n')

fitResults = [rho(:),...
    C11(:),C11FitData(:),...
    C22(:),C22FitData(:),...
    C12(:),C12FitData(:),...
    C33(:),C33FitData(:)];

disp(array2table(fitResults,...
    'VariableNames',{'rho',...
    'C11','C11Fit',...
    'C22','C22Fit',...
    'C12','C12Fit',...
    'C33','C33Fit'}))

fprintf('\n==============================\n')
fprintf('ABSOLUTE FITTING ERRORS\n')
fprintf('==============================\n')

errorC11 = C11FitData(:)-C11(:);
errorC22 = C22FitData(:)-C22(:);
errorC12 = C12FitData(:)-C12(:);
errorC33 = C33FitData(:)-C33(:);

errorResults = [rho(:),errorC11,errorC22,errorC12,errorC33];

disp(array2table(errorResults,...
    'VariableNames',{'rho','errorC11','errorC22','errorC12','errorC33'}))

fprintf('\nMaximum errors:\n')
fprintf('C11 = %.6e\n',max(abs(errorC11)))
fprintf('C22 = %.6e\n',max(abs(errorC22)))
fprintf('C12 = %.6e\n',max(abs(errorC12)))
fprintf('C33 = %.6e\n',max(abs(errorC33)))

dC11Data = df{1,1,1,1}(rho);
dC22Data = df{2,2,2,2}(rho);
dC12Data = df{1,1,2,2}(rho);
dC33Data = df{1,2,1,2}(rho);

fprintf('\n==============================\n')
fprintf('FIRST DERIVATIVES\n')
fprintf('==============================\n')

derivativeResults = [rho(:),...
    dC11Data(:),...
    dC22Data(:),...
    dC12Data(:),...
    dC33Data(:)];

disp(array2table(derivativeResults,...
    'VariableNames',{'rho','dC11','dC22','dC12','dC33'}))

fprintf('\n==============================\n')
fprintf('END POINT CHECK\n')
fprintf('==============================\n')

fprintf('rho = 0\n')
fprintf('C11 = %.12e\n',f{1,1,1,1}(0))
fprintf('C22 = %.12e\n',f{2,2,2,2}(0))
fprintf('C12 = %.12e\n',f{1,1,2,2}(0))
fprintf('C33 = %.12e\n',f{1,2,1,2}(0))

fprintf('\nrho = 1\n')
fprintf('C11 = %.12e\n',f{1,1,1,1}(1))
fprintf('C22 = %.12e\n',f{2,2,2,2}(1))
fprintf('C12 = %.12e\n',f{1,1,2,2}(1))
fprintf('C33 = %.12e\n',f{1,2,1,2}(1))