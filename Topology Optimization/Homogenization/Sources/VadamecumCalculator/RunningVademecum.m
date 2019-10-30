function RunningVademecum

% 
dSmooth = obtainSettings('SmoothRectangleQ2','SmoothRectangle');
dSmooth.smoothingExponentSettings.type = 'Given';
dSmooth.smoothingExponentSettings.q = 2;
computeVademecum(dSmooth);


dSmooth = obtainSettings('RectangleVademecum','Rectangle');
dSmooth.smoothingExponentSettings.type = 'Given';
computeVademecum(dSmooth);

dSmooth = obtainSettings('OptimalSuperEllipseVademecum','SmoothRectangle');
dSmooth.smoothingExponentSettings.type = 'Optimal';
computeVademecum(dSmooth);

end

function computeVademecum(d)
vc = VademecumCellVariablesCalculator(d);
vc.computeVademecumData()
vc.saveVademecumData();
end

function d = obtainSettings(prefix,freeFemFile)
d = SettingsVademecumCellVariablesCalculator();
d.freeFemFileName = freeFemFile;
d.fileName   = prefix;
d.mxMin = 0.01;
d.mxMax = 0.99;
d.myMin = 0.01;
d.myMax = 0.99;
d.nMx   = 20;
d.nMy   = 20;
d.outPutPath = [];
d.print = false;
d.freeFemSettings.hMax = 0.02;%0.0025;
end