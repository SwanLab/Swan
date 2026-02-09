clc,clear,close all
s.monitoring = true;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 250;
s.holeType   = 'ReinforcedHoneycomb';
s.nSteps     = [200];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('HoneycombBenchmark03','mat','phi','degradation')

% save('TestHexagonDensity','mat','phi','holeParam')

% mesh = PFH.InnerMesh;
% mesh.plot
% save('TestHoneycombLevelSet','mat','phi','holeParam','mesh')

s.nu         = 0.3;
s.holeType   = 'Hexagon';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();
[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('HexagonBenchmark03','mat','phi','degradation')

s.nu         = 0.2;
s.holeType   = 'ReinforcedHoneycomb';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();
[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('HoneycombBenchmark02','mat','phi','degradation')

s.nu         = 0.2;
s.holeType   = 'Hexagon';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();
[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('HexagonBenchmark02','mat','phi','degradation')
