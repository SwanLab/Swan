%% RUN
% clc,clear,close all
% s.E          = 1;
% s.nu         = 0.3;
% s.meshType   = 'Square';
% s.meshN      = 150;
% s.holeType   = 'Square';
% s.nSteps     = [2];
% s.pnorm      = 'Inf';
% s.damageType = "Area";
% PFH = TestingPhaseFieldHomogenizer(s);
% [mat,phi,holeParam] = PFH.compute();


%% SAVE + PLOTS 
%save("EllipseMicroDamageArea","mat","phi","holeParam")
figure()
tiledlayout(3,3)
for i=1:9
    nexttile
    i1 = ceil(i/3);
    i2 = rem(i-1,3)+1;
    if ismember(i,[1,2,4,5,9])
    s = surf(holeParam{1},holeParam{2},squeeze(mat(i1,i2,:,:)));
    xlabel('m2')
    ylabel('m1')
    view(2)
    colormap turbo
    %shading interp
    caxis([0,1])
    end
    str = sprintf('C%d%d',i1,i2);
    title(str)
end
sgtitle('Ellipse homogenization')
cb = colorbar;
cb.Layout.Tile = 'East';

% figure()
% s = surf(holeParam{1},holeParam{2},squeeze(phi));
% xlabel('m2')
% ylabel('m1')
% sgtitle('Rectangle Perimeter Damage')