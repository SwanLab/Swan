function RunningOneMicrostructureVademecum


dSmooth = obtainSettings('MySmoothCorner','SmoothRectangle');
vc = VademecumCellVariablesCalculator(dSmooth);
vc.computeVademecumData()

end

function d = obtainSettings(prefix,freeFemFile)
d = SettingsVademecumCellVariablesCalculator();
d.freeFemFileName = freeFemFile;
d.fileName   = prefix;
d.mxMin = 0.5;
d.mxMax = 0.5;
d.myMin = 0.5;
d.myMax = 0.5;
d.nMx   = 1;
d.nMy   = 1;
d.outPutPath = [];
d.print = true;

d.freeFemSettings.hMax  = 0.02;%0.0025;
%d.smoothingExponentSettings.type = 'Optimal';

d.smoothingExponentSettings.type = 'Given';
d.smoothingExponentSettings.q = 2;

end