E = 1;
nu = 0.5;
Gc = 0.1;
l0 =  0.03;

sigC = @(gPrime,cOmega) sqrt(((2*Gc*E)/(l0*cOmega))*((1)/(gPrime)));

gPrimeAT1 = 2;
gPrimeAT2 = 1e20; %infty

lchAT1 = Gc/(0.5*sigC(gPrimeAT1,8/3)^2/E)
lchAT2 = Gc/(0.5*sigC(gPrimeAT2,2)^2/E)

sigCRat = sigC(gPrimeAT1,8/3);
sigCRat = 1;
lchRat = Gc/(0.5*sigCRat^2/E)

fHex = load('HexagonAreaNew.mat');
C11PrimeHexa = fHex.degradation.dfun{1,1,1,1}(0);
C12PrimeHexa = fHex.degradation.dfun{1,1,2,2}(0);
derivFactorHexa = (2*nu*C12PrimeHexa - (1+nu^2)*C11PrimeHexa);
sigCHexa = sqrt(((2*Gc*E)/(l0*(8/3)))*(1/derivFactorHexa));
lchHexa = Gc/(0.5*sigCHexa^2/E)

fHoney = load('HoneycombAreaNew.mat');
C11PrimeHoney = fHoney.degradation.dfun{1,1,1,1}(0);
C12PrimeHoney = fHoney.degradation.dfun{1,1,2,2}(0);
derivFactorHoney = 2*nu*C12PrimeHoney - (1+nu^2)*C11PrimeHoney;
sigCHoney = sqrt(((2*Gc*E)/(l0*(8/3))*(1/(derivFactorHoney))));
lchHoney = Gc/(0.5*sigCHoney^2/E)


k  = E./(2.*(1-nu));
mu = E./(2.*(1+nu));
etak  = mu;
etamu = (k.*mu)./(2.*mu+k);
HSk = (etak/(etak*k+k^2))*k;
HSmu = (etamu/(etamu*mu+mu^2))*mu;


lHSAT1mu = lchAT1*HSmu*(1/(8/3))
lHSAT2mu = lchAT2*HSmu*(1/(2))
lHSRatmu = lchRat*HSmu*(1/(pi))
lHSHexamu = lchHexa*HSmu*(1/(8/3))
lHSHoneymu = lchHoney*HSmu*(1/(8/3))

lHSAT1k = lchAT1*HSk*(1/(8/3))
lHSAT2k = lchAT2*HSk*(1/(2))
lHSRatk = lchRat*HSk*(1/(pi))
lHSHexak = lchHexa*HSk*(1/(8/3))
lHSHoneyk = lchHoney*HSk*(1/(8/3))


