clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Tetradecahedron';
s.meshN      = 2;
s.holeType   = 'Tetradecahedron';%'ReinforcedHoneycomb';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

% [f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
% degradation.fun = f;
% degradation.dfun = df;
% degradation.ddfun = ddf;
save('Tetradecahedron','mat','phi')
