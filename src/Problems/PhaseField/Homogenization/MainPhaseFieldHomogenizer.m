%% Square-perimeter
clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 200;
s.holeType   = 'Square';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = "Perimeter";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
fplot(f,[0 1])
save('SquarePerimeter','mat','phi','degradationFun')

%% Square-area
clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Square';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
fplot(f,[0 1])
save('SquareArea','mat','phi','degradationFun')

%% Circle-perimeter
clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Circle';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = "Perimeter";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
fplot(f,[0 1])
save('CirclePerimeter','mat','phi','degradationFun')

%% Circle-area
clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Circle';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
fplot(f,[0 1])
save('CircleArea','mat','phi','degradationFun')

%% Hexagon-perimeter
clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 100;
s.holeType   = 'SmoothHexagon';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = "Perimeter";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
fplot(f,[0 1])
save('HexagonPerimeter','mat','phi','degradationFun')

%% Hexagon-area
clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 100;
s.holeType   = 'SmoothHexagon';
s.nSteps     = [100];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
degradationFun.fun = f;
degradationFun.dfun = df;
degradationFun.ddfun = ddf;
fplot(f,[0 1])
save('HexagonArea','mat','phi','degradationFun')

%% Ellipse
clc,clear,close all
s.monitoring = false;
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Hexagon';
s.meshN      = 150;
s.holeType   = 'SmoothHexagon';
s.nSteps     = [100,100];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


% [f,df,ddf] = DamageHomogenizationFitter.computePolynomial(9,phi,mat);
% degradationFun.fun = f;
% degradationFun.dfun = df;
% degradationFun.ddfun = ddf;
% fplot(f,[0 1])
save('Ellipse','mat','holeParam')
