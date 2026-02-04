clc,clear,close all
s.monitoring = true;
s.E          = 15;
s.nu         = 0.35;
s.meshType   = 'Hexagon';
s.meshN      = 100;
s.holeType   = 'Hexagon';
s.nSteps     = [1];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizerLevelSet(s);
[mat,phi,holeParam] = PFH.compute();
%save('HexaNu','mat','phi')

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
%save('SmoothHexagonPerle','mat','phi','+degradation')

mesh = PFH.InnerMesh;
mesh.plot
% %save('TestHoneycombLevelSet','mat','phi','holeParam','mesh')
