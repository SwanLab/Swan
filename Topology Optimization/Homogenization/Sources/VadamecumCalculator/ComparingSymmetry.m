function ComparingSymmetry

incPhi = pi/180;


fileName = 'StressSymmetryTraction';
sPnormT = computeExperiment(incPhi,fileName);


fileName = 'StressSymmetryCompression';
sPnormC = computeExperiment(-incPhi,fileName);

end

function sPnorm = computeExperiment(incPhi,fileName)

txi = pi/2 - 0.2;%1083;
rho = 0.9;
q = 4;

phi = 0+incPhi;
pNorm = 'max';
print = false;
hMesh = 0.01;
hasToCaptureImage = true;

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
sN.printStress();



end