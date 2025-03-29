clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 150;
s.holeType   = 'SmoothHexagon';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();

DHF = DamageHomogenizationFitter();
[f,df,ddf] = DHF.computePolynomialFitting(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
save('HexagonMicroDamageArea','mat','phi','degradationFun')