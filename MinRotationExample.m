eps = [1 2; 2 0];
epsV   = [eps(1,1) eps(2,2) 2*eps(1,2)]';
theta = sym("theta");

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

% C = zeros(3,3);
% C(1,1)= 2*mu+l;
% C(1,2)= 10e-10;%l;
% C(2,1)= 10e-10;%l;
% C(2,2)= 10e-10;%2*mu+l;
% C(3,3)= 10e-10;%mu;

Reps = [(1+cos(2*theta))/2 , (1-cos(2*theta))/2 , (-sin(2*theta))/2 ;
        (1-cos(2*theta))/2 , (1+cos(2*theta))/2 , (sin(2*theta))/2  ;
        sin(2*theta)       , -sin(2*theta)      , (cos(2*theta))    ];

Rsig = [(1+cos(2*theta))/2 , (1-cos(2*theta))/2 , sin(2*theta)   ;
        (1-cos(2*theta))/2 , (1+cos(2*theta))/2 , -sin(2*theta)  ;
        (-sin(2*theta))/2  , (sin(2*theta))/2   , (cos(2*theta)) ];

energy = epsV'*Rsig*C*Reps*epsV;

%% Principal directions
[epsDir,epsVal] = eigs(eps);
epsAngle = atan2(epsDir(1,:),epsDir(2,:));
for i=1:length(epsAngle)
    if epsAngle(i)<0
        epsAngle(i) = epsAngle(i)+pi;
    end
end

sigV = C*epsV;
sig = [sigV(1), sigV(3);
       sigV(3), sigV(2)];
[sigDir,sigVal] = eigs(sig);
sigAngle = atan2(sigDir(1,:),sigDir(2,:));
for i=1:length(sigAngle)
    if sigAngle(i)<0
        sigAngle(i) = sigAngle(i)+pi;
    end
end

%% Plots
fplot(theta,energy,[0 pi]);
hold on