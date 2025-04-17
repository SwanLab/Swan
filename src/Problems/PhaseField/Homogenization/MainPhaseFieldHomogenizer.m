clc,clear,close all
s.monitoring = true;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Square';
s.nSteps     = [5];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
fplot(f,[0 1])
%save('HexagonMicroDamageArea','mat','phi','degradationFun')