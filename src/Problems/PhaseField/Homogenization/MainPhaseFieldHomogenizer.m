%% RUN
clc,clear,close all
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 150;
s.holeType   = 'Ellipse';
s.nSteps     = [2,2];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


%% SAVE + PLOTS 
%save("EllipseMicroDamageArea","mat","phi","holeParam")

for i=1:9
    subplot(3,3,i)
    i1 = ceil(i/3);
    i2 = rem(i-1,3)+1;
    if ismember(i,[1,2,4,5,9])
    s = surf(holeParam{1},holeParam{2},squeeze(mat(i1,i2,:,:)));
    xlabel('m2')
    ylabel('m1')
    view(2)
    shading interp
    end
    str = sprintf('C%d%d',i1,i2);
    title(str)
end
sgtitle('Ellipse homogenization')

% figure()
% s = surf(holeParam{1},holeParam{2},squeeze(phi));
% xlabel('m2')
% ylabel('m1')
% sgtitle('Rectangle Perimeter Damage')