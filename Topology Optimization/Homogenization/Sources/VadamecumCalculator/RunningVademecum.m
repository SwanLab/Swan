function RunningVademecum

% 
% dSmooth = obtainSettings('SuperEllipseQMax');
% dSmooth.smoothingExponentSettings.type = 'Given';
% dSmooth.smoothingExponentSettings.q = 32;
% computeVademecum(dSmooth);
% 
% dSmooth = obtainSettings('SuperEllipseQ2');
% dSmooth.smoothingExponentSettings.type = 'Given';
% dSmooth.smoothingExponentSettings.q = 2;
% computeVademecum(dSmooth);

% 
% dSmooth = obtainSettings('RectangleVademecum','Rectangle');
% dSmooth.smoothingExponentSettings.type = 'Given';
% computeVademecum(dSmooth);

dSmooth = obtainSettings('SuperEllipseQOptAnalytic');
dSmooth.smoothingExponentSettings.type = 'Optimal';
computeVademecum(dSmooth);

end

function computeVademecum(d)
vc = VademecumCellVariablesCalculator(d);
vc.computeVademecumData()
vc.saveVademecumData();
end

function d = obtainSettings(prefix)
d = SettingsVademecumCellVariablesCalculator();
d.fileName   = prefix;
d.mxMin = 0.01;
d.mxMax = 0.99;
d.myMin = 0.01;
d.myMax = 0.99;
d.nMx   = 20;
d.nMy   = 20;
d.outPutPath = 'Topology Optimization/Vademecums/';
d.print = false;
end