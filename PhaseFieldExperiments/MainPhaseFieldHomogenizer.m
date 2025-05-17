%% RUN
clc,clear,close all
s.E          = 1;
s.nu         = 0.3;
s.meshType   = 'Square';
s.meshN      = 200;
s.holeType   = 'Rectangle';
s.nSteps     = [200,1];
s.pnorm      = 'Inf';
s.damageType = "Perimeter";
PFH = TestingPhaseFieldHomogenizer(s);
[mat,phi,holeParam] = PFH.compute();
phi = holeParam;
load('CircleMicroDamagePerimeter.mat')
f1 = degradationFun.fun;
load('CirclePerimeter.mat')
for i=1:3
    for j=1:3
        f = degradation.fun{i,j};
        f2{i,j} = @(x) 210.*f(x);
    end
end

%% test
tiledlayout(2,2)
nexttile
hold on
fplot(f1(1,1),[0 1],'Color','#000000','LineWidth',1.5);
fplot(f2(1,1),[0 1],'--','Color','#D95319','LineWidth',1.5);
ylabel(char(8450)+"11 [GPa]");
ylim([0,inf])
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(f1(1,2),[0 1],'Color','#000000','LineWidth',1.5);
fplot(f2(1,2),[0 1],'-','Color','#D95319','LineWidth',1.5);

ylabel(char(8450)+"12 [GPa]");
ylim([0,inf])
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(f1(2,2),[0 1],'Color','#000000','LineWidth',1.5);
fplot(f2(2,2),[0 1],'-','Color','#D95319','LineWidth',1.5);

ylabel(char(8450)+"22 [GPa]");
ylim([0,inf])
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(f1(3,3),[0 1],'Color','#000000','LineWidth',1.5);
fplot(f2(3,3),[0 1],'-','Color','#D95319','LineWidth',1.5);

ylabel(char(8450)+"33 [GPa]");
ylim([0,inf])
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

lg =legend('v1','v2');
lg.Layout.Tile = 'East';