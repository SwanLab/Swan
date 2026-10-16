clc,clear,close all
s.monitoring = true;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 100;
s.holeType   = 'Crack';%'ReinforcedHoneycomb';
s.nSteps     = [3];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

phi = holeParam{1}
%save('HorizontalCrack','mat','phi','holeParam')

% [f,df,ddf] = DamageHomogenizationFitter.computePolynomial(7,holeParam{1},mat);
% degradation.fun = f;
% degradation.dfun = df;
% degradation.ddfun = ddf;
% save('HorizontalCrackCoarse','mat','phi','degradation','holeParam')
