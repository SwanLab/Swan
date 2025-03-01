function StressMeshVariation

hMesh = [0.1, 0.05, 0.0025, 0.00125];
for imesh = 1:length(hMesh)
dSmooth = obtainSettings(['RectangleStressMeshDependency',num2str((imesh))],'Rectangle',hMesh(imesh));
computeVademecum(dSmooth);
end

end

function computeVademecum(d)
vc = VademecumCellVariablesCalculator(d);
vc.computeVademecumData()
end

function d = obtainSettings(prefix,freeFemFile,h)
d = SettingsVademecumCellVariablesCalculator();
d.freeFemFileName = freeFemFile;
d.fileName   = prefix;
d.mxMin = 0.8;
d.mxMax = 0.8;
d.myMin = 0.8;
d.myMax = 0.8;
d.nMx   = 1;
d.nMy   = 1;
d.outPutPath = [];
d.print = true;
d.freeFemSettings.hMax = h;%0.0025;
d.smoothingExponentSettings.type = 'Given';
d.smoothingExponentSettings.q = 2;

end