txi = pi/2 - 0.1083;
rho = 0.9;
q = 4;
phi = 0-pi/6;
pNorm = 'max';
print = true;
hMesh = 0.1;
hasToCaptureImage = true;
fileName = 'StressSymmetry';
mx = SuperEllipseParamsRelator.mx(txi,rho,q); 
my = SuperEllipseParamsRelator.my(txi,rho,q);

s.mx       = mx;
s.my       = my;
s.q        = q;
s.phi      = phi;
s.pNorm    = pNorm;
s.print    = print;
s.hMesh    = hMesh;
s.fileName = fileName;
s.hasToCaptureImage = false;
sN = StressNormSuperEllipseComputer(s);
sPnorm = sN.compute();

