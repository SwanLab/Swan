% clc,clear,close all
% s.E          = 1;
% s.nu         = 0.3;
% s.meshType   = 'Square';
% s.meshN      = 150;
% s.holeType   = 'Ellipse';
% s.nSteps     = [100,100];
% s.pnorm      = 'Inf';
% s.damageType = "Area";
% PFH = PhaseFieldHomogenizer(s);
% [mat,phi,holeParam] = PFH.computeHomogMaterial();
% save("EllipseMicroDamageArea","mat","phi","holeParam")

% clc,clear,close all
% s.E          = 1;
% s.nu         = 0.3;
% s.meshType   = 'Square';
% s.meshN      = 150;
% s.holeType   = 'Rectangle';
% s.nSteps     = [100,100];
% s.pnorm      = 'Inf';
% s.damageType = "Area";
% PFH = PhaseFieldHomogenizer(s);
% [mat,phi,holeParam] = PFH.computeHomogMaterial();
% save("RectangleMicroDamageArea","mat","phi","holeParam")

clc,clear,close all
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Ellipse';
s.nSteps     = [100,100];
s.pnorm      = 'Inf';
s.damageType = "Perimeter";
PFH = PhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.computeHomogMaterial();
save("EllipseMicroDamagePerimeter","mat","phi","holeParam")

clc,clear,close all
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Ellipse';
s.nSteps     = [100,100];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = PhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.computeHomogMaterial();
save("RectangleMicroDamagePerimeter","mat","phi","holeParam")


% for i=1:9
%     subplot(3,3,i)
%     i1 = ceil(i/3);
%     i2 = rem(i-1,3)+1;
%     if ismember(i,[1,2,4,5,9])
%     surf(holeParam{1},holeParam{2},squeeze(mat(i1,i2,:,:)))
%     xlabel('m2')
%     ylabel('m1')
%     end
% end
% title('Ellipse homogenization')