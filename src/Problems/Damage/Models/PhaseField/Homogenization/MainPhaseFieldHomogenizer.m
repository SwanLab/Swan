clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 50;
s.holeType   = 'Ellipse';
s.nSteps     = [2 2];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
fplot(f,[0 1])
%save('CirclePerimeter','mat','phi','degradation')
