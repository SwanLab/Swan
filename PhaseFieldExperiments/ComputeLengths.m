clc,clear
E = 210;
nu = 0.3;
Gc = 2.7e-3;
l0 = 0.025;

C11 = E/((1+nu)*(1-nu));
k  = E./(2.*(1-nu));
mu = E./(2.*(1+nu));

sigC = @(gPrime,cOmega) sqrt(((2*Gc*E)/(l0*cOmega))*((1)/(gPrime)));

gPrimeAT1 = 2;
gPrimeAT2 = 1e20; %infty

lchAT1 = Gc/(0.5*sigC(gPrimeAT1,8/3)^2/E)
lchAT2 = Gc/(0.5*sigC(gPrimeAT2,2)^2/E)

sigCRat = 2.44542;
lchRat = Gc/(0.5*sigCRat^2/E)

fHex = load('HexagonBenchmark03.mat');
C11PrimeHexa = -E*fHex.degradation.dfun{1,1,1,1}(0);
C33PrimeHexa = -E*fHex.degradation.dfun{1,2,1,2}(0);
derivHexaK   = (C11PrimeHexa - C33PrimeHexa)/k;
derivHexaMu  = C33PrimeHexa/mu;
sigCHexaK = sqrt(((2*Gc*k)/(l0*(8/3)))*(1/derivHexaK));
sigCHexaMu = sqrt(((2*Gc*mu)/(l0*(8/3)))*(1/derivHexaMu));
lchHexaK = 2*Gc*k/sigCHexaK^2
lchHexaMu = 2*Gc*mu/sigCHexaMu^2

fHoney = load('HoneycombBenchmark03.mat');
C11PrimeHoney = -E*fHoney.degradation.dfun{1,1,1,1}(0);
C33PrimeHoney = -E*fHoney.degradation.dfun{1,2,1,2}(0);
derivHoneyK   = (C11PrimeHoney - C33PrimeHoney)/k;
derivHoneyMu  = C33PrimeHoney/mu;
sigCHoneyK = sqrt(((2*Gc*k)/(l0*(8/3)))*(1/derivHoneyK));
sigCHoneyMu = sqrt(((2*Gc*mu)/(l0*(8/3)))*(1/derivHoneyMu));
lchHoneyK = 2*Gc*k/sigCHoneyK^2
lchHoneyMu = 2*Gc*mu/sigCHoneyMu^2

etak  = mu;
etamu = (k.*mu)./(2.*mu+k);
HSk = (etak/(etak+k));
HSmu = (etamu/(etamu+mu));


lHSAT1mu = lchAT1*HSmu*(1/(8/3))
lHSAT2mu = lchAT2*HSmu*(1/(2))
lHSRatmu = lchRat*HSmu*(1/(pi))
lHSHexamu = lchHexaMu*HSmu*(1/(8/3))
lHSHoneymu = lchHoneyMu*HSmu*(1/(8/3))

lHSAT1k = lchAT1*HSk*(1/(8/3))
lHSAT2k = lchAT2*HSk*(1/(2))
lHSRatk = lchRat*HSk*(1/(pi))
lHSHexak = lchHexaK*HSk*(1/(8/3))
lHSHoneyk = lchHoneyK*HSk*(1/(8/3))


