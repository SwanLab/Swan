function RunningOneMicrostructureVademecum

%dSmooth = obtainSettings('SquareMesh4b','Rectangle');
dSmooth = obtainSettings('SmallCircleQ2','SmoothRectangle');

vc = VademecumCellVariablesCalculator(dSmooth);
vc.computeVademecumData()

end

function d = obtainSettings(prefix,freeFemFile)
d = SettingsVademecumCellVariablesCalculator();
d.freeFemFileName = freeFemFile;
d.fileName   = prefix;
d.mxMin = 0.1;
d.mxMax = 0.1;
d.myMin = 0.1;
d.myMax = 0.1;
d.nMx   = 1;
d.nMy   = 1;
d.outPutPath = [];
d.print = true;

%i = 4;
%d.freeFemSettings.hMax  = 10^(-1)/(2^(i - 1));%0.0025;

d.freeFemSettings.hMax  = 0.02;%0.0025;

%d.smoothingExponentSettings.type = 'Optimal';

d.smoothingExponentSettings.type = 'Given';
d.smoothingExponentSettings.q = 2;

end