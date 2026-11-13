clc,clear,close all
s.monitoring = true;

s.material.A.E  = 1e-12;
s.material.A.nu = 1e-12;
s.material.B.E  = 1;
s.material.B.nu = 0.3;

s.mesh.type = 'Square';
s.mesh.N    = 250;

s.hole.type = 'CircleCrack';
% s.hole.size = {linspace(4/250,0.999,300)};
s.hole.size = {(4/250):5e-3:0.99};
s.hole.transSize = 0.05;
s.hole.pnorm = 'Inf';


PFH = TestingPhaseFieldHomogenizerV2(s);
[mat,hole] = PFH.compute();

% figure()
% surf(s.hole.size{1},s.hole.size{2},squeeze(mat(1,1,1,1,:,:))')
% xlabel('Crack length (x)')
% ylabel('Crack length (y)')
% 
% figure()
% surf(s.hole.size{1},s.hole.size{2},squeeze(mat(2,2,2,2,:,:))')
% xlabel('Crack length (x)')
% ylabel('Crack length (y)')

phi = hole;
input = s;
save('CircleCrack','mat','phi','input')

[f,df,ddf] = DamageHomogenizationFitter.computePolynomial(7,phi,mat);
degradation.fun = f;
degradation.dfun = df;
degradation.ddfun = ddf;
save('CircleCrack','mat','phi','degradation')
