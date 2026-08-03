clc,clear,close all
%% State (strain)
xi = linspace(0,2*pi,300);
beta = linspace(0,pi,300);
E_eps = zeros(length(xi),length(beta));
BetaOpt = zeros(1,length(xi));
phi = [0 0.01 0.05 0.1];

c2 = @(x) cos(2*x);
s2 = @(x) sin(2*x);
R = @(x) [ (1 + c2(x))/2,  (1 - c2(x))/2,  -s2(x)/sqrt(2);
           (1 - c2(x))/2,  (1 + c2(x))/2,   s2(x)/sqrt(2);
           s2(x)/sqrt(2), -s2(x)/sqrt(2),      c2(x)     ];

eps = @(y) [cos(y); sin(y); 0];

matInfo = load('HorizontalCrack.mat');
funMat = matInfo.degradation.fun;

for i=1:length(phi)
    C = zeros(3,3);
    % Isotropic material
    % mu     = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1,0.3);
    % k      = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1,0.3,2);
    % lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(mu,k,2);
    % C = 2*mu.*eye4D(2) + lambda.*kronEye(2);
    % C = convert2Voigt(C,'Constitutive');

    % Anisotropic material
    C(1,1)= double(subs(funMat{1,1,1,1},phi(i)));
    C(1,2)= double(subs(funMat{1,1,2,2},phi(i)));
    C(2,1)= double(subs(funMat{2,2,1,1},phi(i)));
    C(2,2)= double(subs(funMat{2,2,2,2},phi(i)));
    C(3,3)= double(subs(funMat{1,2,1,2},phi(i)));

    C(3,3) = C(3,3)*2; %Kelvin notation

    for j=1:length(xi)
        for k=1:length(beta)
            xiV   = xi(j);
            betaV = beta(k);
            E_eps(j,k) = 0.5*eps(xiV)'*R(betaV)*C*R(betaV)'*eps(xiV);
        end
        [~,pos] = min(E_eps(j,:));
        BetaOpt(j) = beta(pos);
    end
    figure(2*i-1)
    plot(rad2deg(xi),rad2deg(BetaOpt))
    hold on
    title(['Optimal angle (beta) [phi = ',num2str(phi(i)),']'])
    ylabel('Rotation angle (beta) [deg]')
    ylim([0 180])
    xlabel('Strain ratio angle (xi) [deg]')
    xlim([0 360])

    figure(2*i)
    surf(rad2deg(beta),rad2deg(xi),E_eps)
    title(['Strain energy (Horizontal) [phi = ',num2str(phi(i)),']'])
    ylabel('Strain ratio angle (xi) [deg]')
    ylim([0 360])
    xlabel('Rotation angle (beta) [deg]')
    xlim([0 180])
end


%% State (stress)
clc,clear,close all
xi = linspace(0,2*pi,300);
beta = linspace(0,pi,300);
E_sig = zeros(length(xi),length(beta));
BetaOpt = zeros(1,length(xi));
phi = [0 0.01 0.05 0.1];

c2 = @(x) cos(2*x);
s2 = @(x) sin(2*x);
R = @(x) [ (1 + c2(x))/2,  (1 - c2(x))/2,  -s2(x)/sqrt(2);
           (1 - c2(x))/2,  (1 + c2(x))/2,   s2(x)/sqrt(2);
           s2(x)/sqrt(2), -s2(x)/sqrt(2),      c2(x)     ];

sig = @(y) [cos(y); sin(y); 0];

matInfo = load('HorizontalCrack.mat');
funMat = matInfo.degradation.fun;
dfunMat = matInfo.degradation.dfun;
for i=1:length(phi)
    C = zeros(3,3);
    % Isotropic material
    % mu     = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1,0.3);
    % k      = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1,0.3,2);
    % lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(mu,k,2);
    % C = 2*mu.*eye4D(2) + lambda.*kronEye(2);
    % C = convert2Voigt(C,'Constitutive');

    % Anisotropic material
    C(1,1)= double(subs(funMat{1,1,1,1},phi(i)));
    C(1,2)= double(subs(funMat{1,1,2,2},phi(i)));
    C(2,1)= double(subs(funMat{2,2,1,1},phi(i)));
    C(2,2)= double(subs(funMat{2,2,2,2},phi(i)));
    C(3,3)= double(subs(funMat{1,2,1,2},phi(i)));
    C(3,3) = C(3,3)*2; %Kelvin notation

    dC(1,1)= double(subs(dfunMat{1,1,1,1},phi(i)));
    dC(1,2)= double(subs(dfunMat{1,1,2,2},phi(i)));
    dC(2,1)= double(subs(dfunMat{2,2,1,1},phi(i)));
    dC(2,2)= double(subs(dfunMat{2,2,2,2},phi(i)));
    dC(3,3)= double(subs(dfunMat{1,2,1,2},phi(i)));
    dC(3,3) = dC(3,3)*2; %Kelvin notation

    for j=1:length(xi)
        for k=1:length(beta)
            xiV   = xi(j);
            betaV = beta(k);
            E_sig(j,k) = 0.5*sig(xiV)'*R(betaV)*inv(C)*R(betaV)'*sig(xiV);
            %E_sig(j,k) = 0.5*sig(xiV)'*R(betaV)*inv(C)*dC*inv(C)*R(betaV)'*sig(xiV);
        end
        [~,pos] = max(E_sig(j,:));
        BetaOpt(j) = beta(pos);
    end
    figure(2*i-1)
    plot(rad2deg(xi),rad2deg(BetaOpt))
    hold on
    title(['Optimal angle (beta) [phi = ',num2str(phi(i)),']'])
    ylabel('Rotation angle (beta) [deg]')
    ylim([0 180])
    xlabel('Stress ratio angle (xi) [deg]')
    xlim([0 360])

    figure(2*i)
    surf(rad2deg(beta),rad2deg(xi),E_sig)
    title(['Strain energy derivative (Horizontal) [phi = ',num2str(phi(i)),']'])
    ylabel('Stress ratio angle (xi) [deg]')
    ylim([0 360])
    xlabel('Rotation angle (beta) [deg]')
    xlim([0 180])
end
