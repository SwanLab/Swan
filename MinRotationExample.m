clc,clear,close all
% 
% matInfo = load('HorizontalCrackDamage.mat');
% matInfo.mat = matInfo.mat;
% matInfo.phi  = matInfo.holeParam{1};
% [funMat,dfunMat,ddfunMat] = computeFittingHomogenization(matInfo,1,4);

%% Strain state
% eps = [1 0; 0 0]; % Traction X
% eps = [0 0; 0 1]; % Traction Y
 eps = [0 1; 1 0]; % Pure shear
% eps = [1 0; 0 -1]; % Possion
% eps = [1 0; 0 1]; % Hidrostatic state
% eps = [1 1; 1 1]; % Mixed state1
% eps = [1 2; 2 1]; % Mixed state2
% eps = [2 1; 1 2]; % Mixed state3

epsV  = [eps(1,1) eps(2,2) 2*eps(1,2)]';

%% Rotation matrices
theta = sym("theta");

Reps = [(1+cos(2*theta))/2 , (1-cos(2*theta))/2 , (sin(2*theta))/2  ;
        (1-cos(2*theta))/2 , (1+cos(2*theta))/2 , (-sin(2*theta))/2 ;
        -sin(2*theta)      , sin(2*theta)       , (cos(2*theta))    ];

Rsig = [(1+cos(2*theta))/2 , (1-cos(2*theta))/2 , sin(2*theta)   ;
        (1-cos(2*theta))/2 , (1+cos(2*theta))/2 , -sin(2*theta)  ;
        (-sin(2*theta))/2  , (sin(2*theta))/2   , (cos(2*theta)) ];

%% Material
E  = 1;
nu = 0.3;
k  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(E,nu,2);
mu = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(E,nu);
l  = IsotropicElasticMaterial.computeLambdaFromShearAndBulk(mu,k,2);
C = zeros(3,3);
C(1,1)= 2*mu+l;
C(1,2)= l;
C(2,1)= l;
C(2,2)= 2*mu+l;
C(3,3)= mu;

phi   = linspace(0,1,10);
colors = jet(10);
colors2 = copper(10);
epsAngle = zeros(2,10);
sigAngle = zeros(2,10);

for i=1:10
C = zeros(3,3);
C(1,1)= double(subs(funMat{1,1},phi(i)));
C(1,2)= double(subs(funMat{1,2},phi(i)));
C(2,1)= double(subs(funMat{2,1},phi(i)));
C(2,2)= double(subs(funMat{2,2},phi(i)));
C(3,3)= double(subs(funMat{3,3},phi(i)));

%% Energy values
energy = epsV'*Rsig*C*Reps*epsV;

%% Principal directions
[epsDir,epsVal] = eigs(eps);
epsAngle(:,i) = atan2(epsDir(2,:),epsDir(1,:));
for j=1:size(epsAngle,1)
    if epsAngle(j,i)<0
        epsAngle(j,i) = epsAngle(j,i)+pi;
    end
end

sigV = C*epsV;
sig = [sigV(1), sigV(3);
       sigV(3), sigV(2)];
[sigDir,sigVal] = eigs(sig);
sigAngle(:,i) = atan2(sigDir(2,:),sigDir(1,:));
for k=1:size(sigAngle,1)
    if sigAngle(k,i)<0
        sigAngle(k,i) = sigAngle(k,i)+pi;
    end
end

%% Plots
fplot(theta,energy,[0 pi],'Color',[colors(i,:)],'LineWidth',1.5);
ax = gca;
ax.XTick = 0:pi/8:pi;
ax.XTickLabel = {'0','22.5','45','67.5','90',...
    '112.5','135','157.5','180'};

if i==1
    f = matlabFunction(energy);
    [~,maxV] = fminbnd(@(theta)-f(theta),0,pi);
    ylim([0 -maxV*1.05])
end
%ylim([0,3]) % For hidrostatic case
hold on
%xline([sigAngle(:,i)],'--',["s1-"+i,"s2-"+i],'Color',[colors2(i,:)])
end
xline([sigAngle(:,end)],'--',["s1","s2"])

xl = xline([epsAngle(:,end)],'-',["e1","e2"]);
xl(1).LabelVerticalAlignment   = 'bottom'; xl(2).LabelVerticalAlignment   = 'bottom';

