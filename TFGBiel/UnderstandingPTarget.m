% clear;
% clc;
% close all;

% prob    = TopOptTestTutorialLevelSetNullSpace();
C       = 0.2;
pTarget = 0.2;
alpha   = 0.2;
p       = 1:16;
% desgVar = prob.designVariable;
% mesh    = prob.mesh;

desgVar = rho;

res = zeros(3,size(p,2));

% for ii = 1:size(p,2)
%     s.mesh            = mesh;
%     s.perimeterTarget = pTarget;
%     s.p               = p(ii);
%     s.gradientTest    = LagrangianFunction.create(mesh,1,'P1');
%     perimeter         = PerimeterNormPFunctional(s);
%     [J,dJ]  = perimeter.computeFunctionAndGradient(desgVar);
%     res(1,ii) = (J+1)*pTarget;
% end
% 
% for jj = 1:size(p,2)
%     s.mesh         = mesh;
%     s.alpha        = alpha;
%     s.p            = p(jj);
%     s.gradientTest = LagrangianFunction.create(mesh,1,'P1');
%     isoperimeter    = VolumeNormPFunctional(s);
%     [J,dJ]    = isoperimeter.computeFunctionAndGradient(desgVar);
%     res(2,jj) = (J+1)*alpha;
% end

for jj = 1:size(p,2)
    s.mesh         = mesh;
    s.C            = C;
    s.p            = p(jj);
    s.gradientTest = LagrangianFunction.create(mesh,1,'P1');
    isoperimeter   = IsoPerimetricNormPFunctional(s);
    J         = isoperimeter.computeFunctionAndGradient(desgVar);
    res(3,jj) = (J+1)*C;
end