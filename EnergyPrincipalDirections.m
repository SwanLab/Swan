clc,clear,close all
%% State (strain)
% ratio = linspace(-1,1,50);
% theta = linspace(0,pi,50);
% E_eps = zeros(length(ratio),length(theta));
% phi = [0 0.2 0.4 0.6 0.8 1];
% 
% figure(1)
% tiledlayout(1,length(phi))
% 
% for k=1:length(phi)
% for i=1:length(ratio)
%     for j=1:length(theta)
%     % Define state (strain)
%     r = ratio(i);
%     thetaVal = theta(j);
% 
%      eps = [cos(thetaVal)^2 + r*sin(thetaVal)^2 , (1-r)*sin(thetaVal)*cos(thetaVal);
%            (1-r)*sin(thetaVal)*cos(thetaVal), sin(thetaVal)^2 + r*cos(thetaVal)^2];
% 
%     % eps2 = ((1+r)/2).*eye(2) + ((1-r)/2).*[cos(2*thetaVal)  sin(2*thetaVal);
%     %                                       sin(2*thetaVal) -cos(2*thetaVal)];
% 
%     % Create material
%     % mu     = 1;
%     % lambda = 1;
%     % C = 2*mu.*eye4D(2) + lambda.*kronEye(2);
% 
%     matInfo = load('HorizontalCrack.mat');
%     funMat = matInfo.degradation.dfun;
%     C = zeros(3,3);
%     C(1,1)= double(subs(funMat{1,1,1,1},phi(k)));
%     C(1,2)= double(subs(funMat{1,1,1,2},phi(k)));
%     C(2,1)= double(subs(funMat{2,2,1,1},phi(k)));
%     C(2,2)= double(subs(funMat{2,2,2,2},phi(k)));
%     C(3,3)= double(subs(funMat{1,2,1,2},phi(k)));
%     C = convert2Tensor(C,'Constitutive');
% 
%     % Compute energy
%     sig = tensorprod(C,eps,[3 4],[1 2]);
%     E_eps(i,j) = 0.5*tensorprod(eps,sig,[1 2],[1 2]);
%     end
% end
% 
% nexttile
% title(['Strain energy (Square) [phi = ',num2str(phi(k)),']'])
% surf(theta,ratio,E_eps)
% xlabel('Rotation angle [rad]')
% ylabel('Strain ratio (e2/e1) [-]')
% 
% end



%% State (stress)
ratio = linspace(-1,1,50);
theta = linspace(0,pi,50);
E_sig = zeros(length(ratio),length(theta));
phi = [0 0.2 0.4 0.6 0.8 1];

figure(2)
tiledlayout(1,length(phi))

for k=1:length(phi)
for i=1:length(ratio)
    for j=1:length(theta)
    % Define state (strain)
    r = ratio(i);
    thetaVal = theta(j);

     sig = [cos(thetaVal)^2 + r*sin(thetaVal)^2 , (1-r)*sin(thetaVal)*cos(thetaVal);
           (1-r)*sin(thetaVal)*cos(thetaVal), sin(thetaVal)^2 + r*cos(thetaVal)^2];

    % eps2 = ((1+r)/2).*eye(2) + ((1-r)/2).*[cos(2*thetaVal)  sin(2*thetaVal);
    %                                       sin(2*thetaVal) -cos(2*thetaVal)];
    
    % Create material
    % mu     = 1;
    % lambda = 1;
    % C = 2*mu.*eye4D(2) + lambda.*kronEye(2);

    matInfo = load('HorizontalCrack.mat');
    funMat = matInfo.degradation.dfun;
    C = zeros(3,3);
    C(1,1)= double(subs(funMat{1,1,1,1},phi(k)));
    C(1,2)= double(subs(funMat{1,1,1,2},phi(k)));
    C(2,1)= double(subs(funMat{2,2,1,1},phi(k)));
    C(2,2)= double(subs(funMat{2,2,2,2},phi(k)));
    C(3,3)= double(subs(funMat{1,2,1,2},phi(k)));
    S = convert2Tensor(inv(C),'Compliance');
        
    % Compute energy
    eps = tensorprod(S,sig,[3 4],[1 2]);
    E_sig(i,j) = 0.5*tensorprod(sig,eps,[1 2],[1 2]);
    end
end

nexttile
title(['Stress energy (Horizontal Crack) [phi = ',num2str(phi(k)),']'])
surf(theta,ratio,E_sig)
xlabel('Rotation angle [rad]')
ylabel('Stress ratio (s2/s1) [-]')

end

%% State (stress)
ratio = linspace(-1,1,50);
theta = linspace(0,pi,50);
E_sig = zeros(length(ratio),length(theta));
phi = [0 0.2 0.4 0.6 0.8 1];

figure(3)
tiledlayout(1,length(phi))

for k=1:length(phi)
for i=1:length(ratio)
    for j=1:length(theta)
    % Define state (strain)
    r = ratio(i);
    thetaVal = theta(j);

     sig = [cos(thetaVal)^2 + r*sin(thetaVal)^2 , (1-r)*sin(thetaVal)*cos(thetaVal);
           (1-r)*sin(thetaVal)*cos(thetaVal), sin(thetaVal)^2 + r*cos(thetaVal)^2];

    % eps2 = ((1+r)/2).*eye(2) + ((1-r)/2).*[cos(2*thetaVal)  sin(2*thetaVal);
    %                                       sin(2*thetaVal) -cos(2*thetaVal)];
    
    % Create material
    % mu     = 1;
    % lambda = 1;
    % C = 2*mu.*eye4D(2) + lambda.*kronEye(2);

    matInfo = load('SquareArea.mat');
    funMat = matInfo.degradation.dfun;
    C = zeros(3,3);
    C(1,1)= double(subs(funMat{1,1,1,1},phi(k)));
    C(1,2)= double(subs(funMat{1,1,1,2},phi(k)));
    C(2,1)= double(subs(funMat{2,2,1,1},phi(k)));
    C(2,2)= double(subs(funMat{2,2,2,2},phi(k)));
    C(3,3)= double(subs(funMat{1,2,1,2},phi(k)));
    S = convert2Tensor(inv(C),'Compliance');
        
    % Compute energy
    eps = tensorprod(S,sig,[3 4],[1 2]);
    E_sig(i,j) = 0.5*tensorprod(sig,eps,[1 2],[1 2]);
    end
end

nexttile
title(['Stress energy (Square) [phi = ',num2str(phi(k)),']'])
surf(theta,ratio,E_sig)
xlabel('Rotation angle [rad]')
ylabel('Stress ratio (s2/s1) [-]')

end

%% State (stress)
ratio = linspace(-1,1,50);
theta = linspace(0,pi,50);
E_sig = zeros(length(ratio),length(theta));
phi = [0 0.2 0.4 0.6 0.8 1];

figure(4)
tiledlayout(1,length(phi))

for k=1:length(phi)
for i=1:length(ratio)
    for j=1:length(theta)
    % Define state (strain)
    r = ratio(i);
    thetaVal = theta(j);

     sig = [cos(thetaVal)^2 + r*sin(thetaVal)^2 , (1-r)*sin(thetaVal)*cos(thetaVal);
           (1-r)*sin(thetaVal)*cos(thetaVal), sin(thetaVal)^2 + r*cos(thetaVal)^2];

    % eps2 = ((1+r)/2).*eye(2) + ((1-r)/2).*[cos(2*thetaVal)  sin(2*thetaVal);
    %                                       sin(2*thetaVal) -cos(2*thetaVal)];
    
    % Create material
    % mu     = 1;
    % lambda = 1;
    % C = 2*mu.*eye4D(2) + lambda.*kronEye(2);

    matInfo = load('CircleArea.mat');
    funMat = matInfo.degradation.dfun;
    C = zeros(3,3);
    C(1,1)= double(subs(funMat{1,1,1,1},phi(k)));
    C(1,2)= double(subs(funMat{1,1,1,2},phi(k)));
    C(2,1)= double(subs(funMat{2,2,1,1},phi(k)));
    C(2,2)= double(subs(funMat{2,2,2,2},phi(k)));
    C(3,3)= double(subs(funMat{1,2,1,2},phi(k)));
    S = convert2Tensor(inv(C),'Compliance');
        
    % Compute energy
    eps = tensorprod(S,sig,[3 4],[1 2]);
    E_sig(i,j) = 0.5*tensorprod(sig,eps,[1 2],[1 2]);
    end
end

nexttile
title(['Stress energy (Circle) [phi = ',num2str(phi(k)),']'])
surf(theta,ratio,E_sig)
xlabel('Rotation angle [rad]')
ylabel('Stress ratio (s2/s1) [-]')

end