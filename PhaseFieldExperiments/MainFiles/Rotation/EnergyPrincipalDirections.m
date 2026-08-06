clc,clear,close all
%% State (strain)
xi = linspace(0,2*pi,17);
beta = linspace(0,pi,300);
phi = [0 0.01 0.05 0.1 0.3];

E_eps = zeros(length(xi),length(beta),length(phi));
EOpt_eps = zeros(length(xi),length(phi));
BetaOptEps = zeros(length(xi),length(phi));
theta_eps = zeros(length(xi),length(phi));
theta_sig = zeros(length(xi),length(phi));
diffAngle = zeros(length(xi),length(phi));


cmap = hsv(length(xi));

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
    mu     = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1,0.3);
    k      = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1,0.3,2);
    lambda = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(mu,k,2);
    Ciso = 2*mu.*eye4D(2) + lambda.*kronEye(2);
    Ciso = convert2Voigt(Ciso,'Constitutive');
    Ciso(3,3) = Ciso(3,3)*2; %Kelvin notation

    % Anisotropic material
    C(1,1)= double(subs(funMat{1,1,1,1},phi(i)));
    C(1,2)= double(subs(funMat{1,1,2,2},phi(i)));
    C(2,1)= double(subs(funMat{2,2,1,1},phi(i)));
    C(2,2)= double(subs(funMat{2,2,2,2},phi(i)));
    C(3,3)= double(subs(funMat{1,2,1,2},phi(i)));

    C(3,3) = C(3,3)*2; %Kelvin notation

    for j=1:length(xi)
        xiV  = xi(j);
        epsV = eps(xiV);
        for k=1:length(beta)
            betaV = beta(k);
            E_eps(j,k,i) = 0.5*epsV'*R(betaV)*C*R(betaV)'*epsV;
        end
        idx = find(abs(E_eps(j,:,i)-min(E_eps(j,:,i)))<1e-10);
        BetaOptEps(j,i) = min(beta(idx));

        sigV = R(BetaOptEps(j,i))*C*R(BetaOptEps(j,i))'*epsV;
        %sigV = Ciso*epsV;
        theta_eps(j,i) = 0.5*atan2d(2*epsV(3),epsV(1)-epsV(2));
        theta_sig(j,i) = 0.5*atan2d(2*sigV(3)/sqrt(2),sigV(1)-sigV(2));

        a_eps = epsV(1) - epsV(2); a_sig = sigV(1) - sigV(2);
        b_eps = 2*epsV(3); b_sig = 2*sigV(3);
        diffAngle(j,i) = 0.5*atan2d(b_sig*a_eps - a_sig*b_eps, a_sig*a_eps + b_sig*b_eps);

        EOpt_eps(j,i) = E_eps(j,min(idx),i);
    end
    

    % figure(3*i-2)
    % plot(rad2deg(xi),theta_eps)
    % title(['Angle principal directions (strain) (Horizontal) [phi = ',num2str(phi(i)),']'])
    % ylabel('Angle principal directions of strain (theta eps) [deg]')
    % ylim([0 360])
    % xlabel('Strain ratio angle (xi) [deg]')
    % xlim([0 360])

    % figure(3*i-1)
    % plot(rad2deg(xi),theta_sig)
    % title(['Angle between principal directions (stress) (Horizontal) [phi = ',num2str(phi(i)),']'])
    % ylabel('Angle principal directions of stress (theta sig) [deg]')
    % ylim([0 360])
    % xlabel('Strain ratio angle (xi) [deg]')
    % xlim([0 360])
    
    % figure(3*i)
    % plot(rad2deg(xi),diffAngle)
    % title(['Angle between principal directions (stress/strain) (Horizontal) [phi = ',num2str(phi(i)),']'])
    % ylabel('Angle between principal directions (Delta theta) [deg]')
    % ylim([0 360])
    % xlabel('Strain ratio angle (xi) [deg]')
    % xlim([0 360])
    
    %Plots
    % figure(4*i-3)
    % plot(rad2deg(xi),rad2deg(BetaOpt))
    % hold on
    % title(['Optimal angle (beta) [phi = ',num2str(phi(i)),']'])
    % ylabel('Rotation angle (beta) [deg]')
    % ylim([0 180])
    % xlabel('Strain ratio angle (xi) [deg]')
    % xlim([0 360])
    % 
    % figure(4*i-2)
    % surf(rad2deg(beta),rad2deg(xi),E_eps)
    % title(['Strain energy (Horizontal) [phi = ',num2str(phi(i)),']'])
    % ylabel('Strain ratio angle (xi) [deg]')
    % ylim([0 360])
    % xlabel('Rotation angle (beta) [deg]')
    % xlim([0 180])
    % 
    % figure(4*i-1)
    % grid on
    % hold on
    % for j=1:length(xi)
    %     plot3(repmat(rad2deg(xi(j)),1,length(beta)),rad2deg(beta),E_eps(j,:),'Color',cmap(j,:))
    % end
    % % lgd = legend('0','22.5','45','67.5','90','112.5','135','157.5','180','202.5','225','247.5','270','292.5','315','337.5','360');
    % % lgd.AutoUpdate = 'off'; 
    % for j=1:length(xi)
    %     idx = find(abs(E_eps(j,:)-min(E_eps(j,:)))<1e-10);
    %     scatter3(repmat(rad2deg(xi(j)),1,length(idx)),rad2deg(beta(idx)),E_eps(j,idx),'o','MarkerEdgeColor',cmap(j,:))
    % end
    % title(['Strain energy (Horizontal) [phi = ',num2str(phi(i)),']'])
    % xlabel('Strain ratio angle (xi) [deg]')
    % xlim([0 360])
    % ylabel('Rotation angle (beta) [deg]')
    % ylim([0 180])
    % view(348,6)
    % 
    % figure(4*i)
    % grid on
    % hold on
    % for j=1:length(xi)
    %     idx = find(abs(E_eps(j,:)-min(E_eps(j,:)))<1e-10);
    %     scatter(repmat(rad2deg(xi(j)),1,length(idx)),rad2deg(beta(idx)),'MarkerEdgeColor','k')
    % end   
    % title(['Optimal angle (beta) [phi = ',num2str(phi(i)),']'])
    % ylabel('Rotation angle (beta) [deg]')
    % ylim([0 180])
    % xlabel('Strain ratio angle (xi) [deg]')
    % xlim([0 360])
end



%% State (stress)
E_sig = zeros(length(xi),length(beta),length(phi));
EOpt_sig = zeros(length(xi),length(phi));
BetaOptSig = zeros(length(xi),length(phi));

cmap = hsv(length(xi));

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
            E_sig(j,k,i) = 0.5*sig(xiV)'*R(betaV)*inv(C)*R(betaV)'*sig(xiV);
            %E_sig(j,k) = 0.5*sig(xiV)'*R(betaV)*inv(C)*dC*inv(C)*R(betaV)'*sig(xiV);
        end
        idx = find(abs(E_sig(j,:,i)-max(E_sig(j,:,i)))<1e-10);
        BetaOptSig(j,i) = min(beta(idx));

        EOpt_sig(j,i) = E_eps(j,min(idx),i);
    end

    figure(4*i-3)
    hold on
    plot(rad2deg(xi),rad2deg(BetaOptSig(:,i)))
    plot(rad2deg(xi),rad2deg(BetaOptEps(:,i)))
    title(['Optimal angle (beta) [phi = ',num2str(phi(i)),']'])
    ylabel('Rotation angle (beta) [deg]')
    ylim([0 180])
    xlabel('Stress ratio angle (xi) [deg]')
    xlim([0 360])
    legend('Beta stress','Beta strain')
    % 
    % figure(4*i-2)
    % surf(rad2deg(beta),rad2deg(xi),E_sig(:,:,i))
    % title(['Stress energy derivative (Horizontal) [phi = ',num2str(phi(i)),']'])
    % ylabel('Stress ratio angle (xi) [deg]')
    % ylim([0 360])
    % xlabel('Rotation angle (beta) [deg]')
    % xlim([0 180])
    % 
    % figure(4*i-1)
    % grid on
    % hold on
    % for j=1:length(xi)
    %     plot3(repmat(rad2deg(xi(j)),1,length(beta)),rad2deg(beta),E_sig(j,:,i),'Color',cmap(j,:))
    % end
    % % lgd = legend('0','22.5','45','67.5','90','112.5','135','157.5','180','202.5','225','247.5','270','292.5','315','337.5','360');
    % % lgd.AutoUpdate = 'off'; 
    % for j=1:length(xi)
    %     idx = find(abs(E_sig(j,:,i)-max(E_sig(j,:,i)))<1e-10);
    %     scatter3(repmat(rad2deg(xi(j)),1,length(idx)),rad2deg(beta(idx)),E_sig(j,idx,i),'o','MarkerEdgeColor',cmap(j,:))
    % end
    % title(['Stress energy (Horizontal) [phi = ',num2str(phi(i)),']'])
    % xlabel('Stress ratio angle (xi) [deg]')
    % xlim([0 360])
    % ylabel('Rotation angle (beta) [deg]')
    % ylim([0 180])
    % view(348,6)
    % 
    % figure(4*i)
    % grid on
    % hold on
    % for j=1:length(xi)
    %     idx = find(abs(E_sig(j,:,i)-max(E_sig(j,:,i)))<1e-10);
    %     scatter(repmat(rad2deg(xi(j)),1,length(idx)),rad2deg(beta(idx)),'MarkerEdgeColor','k')
    % end   
    % title(['Optimal angle (beta) [phi = ',num2str(phi(i)),']'])
    % ylabel('Rotation angle (beta) [deg]')
    % ylim([0 180])
    % xlabel('Stress ratio angle (xi) [deg]')
    % xlim([0 360])
    
end



    figure(100)
    title(['Strain energy at optimal beta [phi = ',num2str(phi(5)),']'])
    hold on
    plot(rad2deg(xi),EOpt_eps(:,5))
    plot(rad2deg(xi),EOpt_sig(:,5),'--')
    xlabel('Ratio angle (xi) [deg]')
    xlim([0 360])
    legend('Beta strain optimal', 'Beta stress optimal')

    abs(EOpt_eps(:,5)-EOpt_sig(:,5))