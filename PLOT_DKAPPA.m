nu = linspace(-1,1,100);
matNu = zeros(2,2,2,2,2,length(nu));
for i=1:length(nu)
    s.monitoring = false;
    s.E          = 210;
    s.nu         = nu(i);
    s.meshType   = 'Hexagon';
    s.meshN      = 100;
    s.holeType   = 'Hexagon';
    s.nSteps     = [2];
    s.pnorm      = 'Inf';
    s.damageType = 'Area';
    PFH = TestingPhaseFieldHomogenizer(s);
    [mat,phi,holeParam] = PFH.compute();
    matNu(:,:,:,:,:,i) = mat(:,:,:,:,:);
end

dmuHomog = zeros(1,length(nu));
dkHomog = zeros(1,length(nu));
for j=1:length(nu)
    mu0 = (1/2).*matNu(1,2,1,2,1,j);
    mu1 = (1/2).*matNu(1,2,1,2,2,j);
    k0  = matNu(1,1,1,1,1,j) - mu0;
    k1  = matNu(1,1,1,1,2,j) - mu1;
    dmuHomog(j) = -(mu0-mu1)/phi(2);
    dkHomog(j)  = -(k0-k1)/phi(2);
end



%% Functions

load('homogNu.mat')
nu = linspace(-1,1,100);


dmuAT1 = zeros(1,length(nu));
dkAT1  = zeros(1,length(nu));
dmuAT2 = zeros(1,length(nu));
dkAT2  = zeros(1,length(nu));
dkRatio = zeros(1,length(nu));
dmuRatio = zeros(1,length(nu));
dkGamma = zeros(1,length(nu));
dmuGamma = zeros(1,length(nu));
dkUB = zeros(1,length(nu));
dkSIMPALL = zeros(1,length(nu));

syms phi
for x=1:length(nu)
E  = 210;
k  = E./(2.*(1-nu(x)));
mu = E./(2.*(1+nu(x)));

dkAT1(x)  = -2*k;
dmuAT1(x) = -2*mu;

k0     = 1e-10;
k1     = k;
mu0    = 1e-10;
mu1    = mu;
etak0  = mu0;
etak1  = mu1;
etamu0 = (k0.*mu0)./(2.*mu0+k0);
etamu1 = (k1.*mu1)./(2.*mu1+k1);

muSIMPALL = ((mu1-mu0).*(etamu0-etamu1).*(phi).*(1-phi) + mu0.*(mu1+etamu0).*(phi) + mu1.*(mu0+etamu1).*(1-phi))./...
                      ((mu1+etamu0).*(phi) + (mu0+etamu1).*(1-phi));

kSIMPALL = ((k1-k0).*(etak0-etak1).*(phi).*(1-phi) + k0.*(k1+etak0).*(phi) + k1.*(k0+etak1).*(1-phi))./...
                         ((k1+etak0).*(phi) + (k0+etak1).*(1-phi));
dkSIMPALL(x) = double(subs(diff(kSIMPALL),0));

kUB  = k0.*(phi) + k1.*(1-phi) - ((1-phi).*phi.*(k1-k0).^2)./(k0.*(1-phi) + k1.*phi + etak1);
dkUB(x) = double(subs(diff(kUB),0));


% Rational (Wu)
a1 = -6;
kRatio  = k.*(phi.^2 - 2.*phi + 1)./(1-(a1+2).*phi);
muRatio = mu.*(phi.^2 - 2.*phi + 1)./(1-(a1+2).*phi);

dkRatio(x) = double(subs(diff(kRatio),0));
dmuRatio(x) = double(subs(diff(muRatio),0));

% Rational (Alessi)
gammaK  = 5;
gammaMu = 1.5;
w = 1 - (1-phi).^2;
kGamma = k.*(1-w)/(1+(gammaK-1)*w);
muGamma = mu.*(1-w)/(1+(gammaMu-1)*w);

dkGamma(x) = double(subs(diff(kGamma),0));
dmuGamma(x) = double(subs(diff(muGamma),0));

end

%% Plots
close all
figure(1)
hold on
plot(nu(2:end-1),dkAT1(2:end-1))
plot(nu(2:end-1),dkRatio(2:end-1))
plot(nu(2:end-1),dkGamma(2:end-1))
plot(nu(2:end-1),dkSIMPALL(2:end-1))
plot(nu(2:end-1),dkUB(2:end-1))
%plot(nu(2:end-1),dkHomog(2:end-1))
legend('AT1','Wu','Alessi','SIMPALL','UB','Homog')


