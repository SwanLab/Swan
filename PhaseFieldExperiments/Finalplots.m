%% Constitutive tensor (after running C_Fun_Ploy)
tiledlayout(1,3)
nexttile
hold on
fplot(funMat(1,1,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(funMat(1,1,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(funMat(1,1,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(funMat(1,1,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(funMat(1,1,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8450)+"11 [GPa]");
ylim([0,inf])
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(funMat(1,2,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(funMat(1,2,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(funMat(1,2,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(funMat(1,2,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(funMat(1,2,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8450)+"12 [GPa]");
ylim([0,inf])
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(funMat(3,3,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(funMat(3,3,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(funMat(3,3,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(funMat(3,3,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(funMat(3,3,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8450)+"33 [GPa]");
ylim([0,inf])
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

lg =legend('Analytical','Circle (Volume)','Circle (Length)','Square (Volume)','Square (Length)');
lg.Layout.Tile = 'East';

%% Constitutive tensor derivative (after running C_Fun_Ploy)
tiledlayout(1,3)
nexttile
hold on
fplot(dfunMat(1,1,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(dfunMat(1,1,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(dfunMat(1,1,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(dfunMat(1,1,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(dfunMat(1,1,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8706)+""+char(8450)+"11/"+char(8706)+char(632)+" [GPa]");
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(dfunMat(1,2,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(dfunMat(1,2,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(dfunMat(1,2,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(dfunMat(1,2,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(dfunMat(1,2,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8706)+""+char(8450)+"12/"+char(8706)+char(632)+" [GPa]");
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(dfunMat(3,3,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(dfunMat(3,3,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(dfunMat(3,3,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(dfunMat(3,3,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(dfunMat(3,3,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8706)+""+char(8450)+"33/"+char(8706)+char(632)+" [GPa]");
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

lg =legend('Analytical','Circle (Volume)','Circle (Length)','Square (Volume)','Square (Length)');
lg.Layout.Tile = 'East';

%% Constitutive tensor second derivative (after running C_Fun_Ploy)
tiledlayout(1,3)
nexttile
hold on
fplot(ddfunMat(1,1,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(ddfunMat(1,1,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(ddfunMat(1,1,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(ddfunMat(1,1,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(ddfunMat(1,1,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8706)+""+char(178)+char(8450)+"11/"+char(8706)+char(632)+char(178)+" [GPa]");
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(ddfunMat(1,2,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(ddfunMat(1,2,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(ddfunMat(1,2,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(ddfunMat(1,2,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(ddfunMat(1,2,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8706)+""+char(178)+char(8450)+"12/"+char(8706)+char(632)+char(178)+" [GPa]");
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

nexttile
hold on
fplot(ddfunMat(3,3,5),[0 1],'Color','#000000','LineWidth',1.5);
fplot(ddfunMat(3,3,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
fplot(ddfunMat(3,3,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
fplot(ddfunMat(3,3,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
fplot(ddfunMat(3,3,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
ylabel(char(8706)+""+char(178)+char(8450)+"33/"+char(8706)+char(632)+char(178)+" [GPa]");
xlabel("Damage "+char(632)+" [-]");
fontsize(gcf,25,'points')

lg =legend('Analytical','Circle (Volume)','Circle (Length)','Square (Volume)','Square (Length)');
lg.Layout.Tile = 'East';

%% Force displacement plots (after loading all outputs)

%% Damage displacement plots

%% Damage fields