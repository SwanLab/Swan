clc,clear,close all

H = TutorialHomogenizationHoleSize();

rho = H.rho;
mat = H.Chomog;

degree = 8;

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(degree,rho,mat);

homogenization.fun = f;
homogenization.dfun = df;
homogenization.ddfun = ddf;

save(fullfile('Gabriel','HMVademecum','Homogenization','SquareHoleSize'),...
    'mat','rho','homogenization')
rhoPlot = linspace(0,1,200);

C11Fit = f{1,1,1,1}(rhoPlot);
C22Fit = f{2,2,2,2}(rhoPlot);
C12Fit = f{1,1,2,2}(rhoPlot);
C33Fit = f{1,2,1,2}(rhoPlot);

dC11Fit = df{1,1,1,1}(rhoPlot);
dC22Fit = df{2,2,2,2}(rhoPlot);
dC12Fit = df{1,1,2,2}(rhoPlot);
dC33Fit = df{1,2,1,2}(rhoPlot);

figure
plot(rhoPlot,C11Fit)
xlabel('\rho_h')
ylabel('C_{1111}^h')
grid on

figure
plot(rhoPlot,C22Fit)
xlabel('\rho_h')
ylabel('C_{2222}^h')
grid on

figure
plot(rhoPlot,C12Fit)
xlabel('\rho_h')
ylabel('C_{1122}^h')
grid on

figure
plot(rhoPlot,C33Fit)
xlabel('\rho_h')
ylabel('C_{1212}^h')
grid on

figure
plot(rhoPlot,dC11Fit)
xlabel('\rho_h')
ylabel('\partial C_{1111}^h/\partial \rho_h')
grid on

figure
plot(rhoPlot,dC33Fit)
xlabel('\rho_h')
ylabel('\partial C_{1212}^h/\partial \rho_h')
grid on