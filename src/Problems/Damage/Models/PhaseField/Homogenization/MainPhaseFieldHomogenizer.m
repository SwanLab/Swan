clc,clear,close all
s.monitoring = true;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 250;
s.holeType   = 'Hexagon';%'ReinforcedHoneycomb';
s.nSteps     = [200];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('HexagonDerivativeNu','mat','phi')
