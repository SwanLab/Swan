firstPart  = fullfile( '/home','alex','Desktop');
secondPart = fullfile('TheFolder/');
outPutPath = fullfile(firstPart,secondPart);
d.fileName = 'VademecumSmoothCorner';
d.outPutPath = fullfile(outPutPath,[d.fileName,'/']);
vc = VademecumComputerAndPlotter(d);
vc.compute();
