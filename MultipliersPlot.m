
%% 4x4 Density

clear;
clc;

%
p0Mesh = QuadMesh(1,1,4,4);
fig = openfig('Paper/Domain4x4/MonitoringDensityGripper4x4Segment1.1.fig');
%

l1 = mean(fig.Children(11).Children.YData(end-50:end));
l2 = mean(fig.Children(15).Children.YData(end-50:end));
l3 = mean(fig.Children(19).Children.YData(end-50:end));
l4 = mean(fig.Children(23).Children.YData(end-50:end));

l5 = mean(fig.Children(10).Children.YData(end-50:end));
l6 = mean(fig.Children(14).Children.YData(end-50:end));
l7 = mean(fig.Children(18).Children.YData(end-50:end));
l8 = mean(fig.Children(22).Children.YData(end-50:end));

l9  = mean(fig.Children(9).Children.YData(end-50:end));
l10 = mean(fig.Children(13).Children.YData(end-50:end));
l11 = mean(fig.Children(17).Children.YData(end-50:end));
l12 = mean(fig.Children(21).Children.YData(end-50:end));

l13 = mean(fig.Children(8).Children.YData(end-50:end));
l14 = mean(fig.Children(12).Children.YData(end-50:end));
l15 = mean(fig.Children(16).Children.YData(end-50:end));
l16 = mean(fig.Children(20).Children.YData(end-50:end));

close all

fL = LagrangianFunction.create(p0Mesh,1,'P0');
fL.setFValues([l1;l2;l3;l4;l5;l6;l7;l8;l9;l10;l11;l12;l13;l14;l15;l16]);
fL.print('Prova');



%% 4x4 Level-set

clear;
clc;

%
p0Mesh = QuadMesh(1,1,4,4);
fig = openfig('Paper/Domain4x4/MonitoringLevelSetGripper4x4Segment1.1.fig');
%

l1 = mean(fig.Children(14).Children.YData(end-50:end));
l2 = mean(fig.Children(18).Children.YData(end-50:end));
l3 = mean(fig.Children(22).Children.YData(end-50:end));
l4 = mean(fig.Children(26).Children.YData(end-50:end));

l5 = mean(fig.Children(13).Children.YData(end-50:end));
l6 = mean(fig.Children(17).Children.YData(end-50:end));
l7 = mean(fig.Children(21).Children.YData(end-50:end));
l8 = mean(fig.Children(25).Children.YData(end-50:end));

l9  = mean(fig.Children(12).Children.YData(end-50:end));
l10 = mean(fig.Children(16).Children.YData(end-50:end));
l11 = mean(fig.Children(20).Children.YData(end-50:end));
l12 = mean(fig.Children(24).Children.YData(end-50:end));

l13 = mean(fig.Children(11).Children.YData(end-50:end));
l14 = mean(fig.Children(15).Children.YData(end-50:end));
l15 = mean(fig.Children(19).Children.YData(end-50:end));
l16 = mean(fig.Children(23).Children.YData(end-50:end));

close all

fL = LagrangianFunction.create(p0Mesh,1,'P0');
fL.setFValues([l1;l2;l3;l4;l5;l6;l7;l8;l9;l10;l11;l12;l13;l14;l15;l16]);
fL.print('Prova');