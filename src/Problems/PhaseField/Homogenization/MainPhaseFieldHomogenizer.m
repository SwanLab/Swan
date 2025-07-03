clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 200;
s.holeType   = 'Circle';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = 'Area';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
fplot(f,[0 1])
save('CircleArea4T','mat','phi','degradation')

s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 100;
s.holeType   = 'Circle';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = 'Perimeter';
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
fplot(f,[0 1])
save('CirclePerimeter4T','mat','phi','degradation')
