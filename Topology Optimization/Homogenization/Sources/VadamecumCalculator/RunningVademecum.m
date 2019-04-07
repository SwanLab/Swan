prefixName = 'VademecumSmoothCorner';
d = SettingsVademecumCellVariablesCalculator();
d.fileName   = prefixName;
d.freeFemFileName = 'SmoothRectangle';
d.mxMin = 0.01;
d.mxMax = 0.99;
d.myMin = 0.01;
d.myMax = 0.99;
d.nMx   = 20;
d.nMy   = 20;
d.outPutPath = [];
d.print = true;
d.freeFemSettings.hMax = 0.02;%0.0025;
vc = VademecumCellVariablesCalculator(d);
vc.computeVademecumData()
vc.saveVademecumData();


% prefixName = 'VademecumCorner';
% d = SettingsVademecumCellVariablesCalculator();
% d.fileName   = prefixName;
% d.freeFemFileName = 'Rectangle';
% d.mxMin = 0.01;
% d.mxMax = 0.99;
% d.myMin = 0.01;
% d.myMax = 0.99;
% d.nMx   = 20;
% d.nMy   = 20;
% d.outPutPath = [];
% d.print = true;
% d.freeFemSettings.hMax = 0.02;%0.0025;
% vc = VademecumCellVariablesCalculator(d);
% vc.computeVademecumData()
% vc.saveVademecumData();