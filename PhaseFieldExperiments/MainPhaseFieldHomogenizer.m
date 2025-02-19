%% RUN
clc,clear,close all
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 200;
s.holeType   = 'Rectangle';
s.nSteps     = [10,10];
s.pnorm      = 'Inf';
s.damageType = "Area";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();


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
    %view([0,0,180])
    view([45,45,45])
    colormap turbo
    caxis([0,1])
    end
    str = sprintf('C%d%d',i1,i2);
    title(str)
end
sgtitle('Ellipse homogenization')
cb = colorbar;
cb.Layout.Tile = 'East';

% surf(holeParam{1},holeParam{2},squeeze(mat(1,1,:,:)));
% xlabel('m2')
% ylabel('m1')
% view(2)
% colormap turbo
% caxis([0,1])
% title(char(8450)+"11")
% cb = colorbar;
% fontsize(gcf,30,'points')

% figure()
% s = surf(holeParam{1},holeParam{2},squeeze(phi));
% xlabel('m2')
% ylabel('m1')
% sgtitle('Rectangle Perimeter Damage')