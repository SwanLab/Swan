%% Comparative constitutive tensor (data vs. fitting)
% close all
% matIdx = 6;
% 
% tiledlayout(2,2)
% nexttile
% hold on
% scatter(matType{matIdx}.phi,squeeze(matType{matIdx}.mat(1,1,:)),'X','LineWidth',1.5)
% fplot(funMat(1,1,matIdx),[0 1],'LineWidth',1.5)
% ylabel(char(8450)+"11 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% %fontsize(gcf,25,'points')
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% 
% nexttile
% hold on
% scatter(matType{matIdx}.phi,squeeze(matType{matIdx}.mat(1,2,:)),'X','LineWidth',1.5)
% fplot(funMat(1,2,matIdx),[0 1],'LineWidth',1.5)
% ylabel(char(8450)+"12 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% %fontsize(gcf,25,'points')
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% 
% nexttile
% hold on
% scatter(matType{matIdx}.phi,squeeze(matType{matIdx}.mat(2,2,:)),'X','LineWidth',1.5)
% fplot(funMat(2,2,matIdx),[0 1],'LineWidth',1.5)
% ylabel(char(8450)+"22 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% %fontsize(gcf,25,'points')
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% 
% nexttile
% hold on
% scatter(matType{matIdx}.phi,squeeze(matType{matIdx}.mat(3,3,:)),'X','LineWidth',1.5)
% fplot(funMat(3,3,matIdx),[0 1],'LineWidth',1.5)
% ylabel(char(8450)+"33 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% %fontsize(gcf,25,'points')
% % set(gca,'xscale','log')
% % set(gca,'yscale','log')
% 
% lg =legend('Data','Fitting');
% lg.Layout.Tile = 'East';

% %% Constitutive tensor (after running C_Fun_Ploy)
% close all
% tiledlayout(2,2)
% nexttile
% hold on
% fplot(funMat(1,1,5),[0 1],'-','Color','#000000','LineWidth',1.5);
% fplot(funMat(1,1,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(funMat(1,1,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(funMat(1,1,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(1,1,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(1,1,6),[0 1],'--','Color',"#7E2F8E",'LineWidth',1.5);
% ylabel(char(8450)+"11 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% nexttile
% hold on
% fplot(funMat(1,2,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(funMat(1,2,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(funMat(1,2,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(funMat(1,2,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(1,2,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(1,2,6),[0 1],'--','Color',"#7E2F8E",'LineWidth',1.5);
% ylabel(char(8450)+"12 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% nexttile
% hold on
% fplot(funMat(2,2,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(funMat(2,2,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(funMat(2,2,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(funMat(2,2,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(2,2,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(2,2,6),[0 1],'--','Color',"#7E2F8E",'LineWidth',1.5);
% ylabel(char(8450)+"22 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% nexttile
% hold on
% fplot(funMat(3,3,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(funMat(3,3,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(funMat(3,3,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(funMat(3,3,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(3,3,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% fplot(funMat(3,3,6),[0 1],'--','Color',"#7E2F8E",'LineWidth',1.5);
% ylabel(char(8450)+"33 [GPa]");
% ylim([0,inf])
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% lg =legend('Analytical','Circle (Volume)','Circle (Length)','Square (Volume)','Square (Length)');
% lg.Layout.Tile = 'East';
% 
% %% Constitutive tensor derivative (after running C_Fun_Ploy)
% tiledlayout(1,3)
% nexttile
% hold on
% fplot(dfunMat(1,1,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(dfunMat(1,1,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(dfunMat(1,1,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(dfunMat(1,1,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(dfunMat(1,1,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% ylabel(char(8706)+""+char(8450)+"11/"+char(8706)+char(632)+" [GPa]");
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% nexttile
% hold on
% fplot(dfunMat(1,2,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(dfunMat(1,2,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(dfunMat(1,2,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(dfunMat(1,2,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(dfunMat(1,2,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% ylabel(char(8706)+""+char(8450)+"12/"+char(8706)+char(632)+" [GPa]");
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% nexttile
% hold on
% fplot(dfunMat(3,3,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(dfunMat(3,3,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(dfunMat(3,3,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(dfunMat(3,3,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(dfunMat(3,3,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% ylabel(char(8706)+""+char(8450)+"33/"+char(8706)+char(632)+" [GPa]");
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% lg =legend('Analytical','Circle (Volume)','Circle (Length)','Square (Volume)','Square (Length)');
% lg.Layout.Tile = 'East';
% 
% %% Constitutive tensor second derivative (after running C_Fun_Ploy)
% tiledlayout(1,3)
% nexttile
% hold on
% fplot(ddfunMat(1,1,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(ddfunMat(1,1,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(ddfunMat(1,1,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(ddfunMat(1,1,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(ddfunMat(1,1,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% ylabel(char(8706)+""+char(178)+char(8450)+"11/"+char(8706)+char(632)+char(178)+" [GPa]");
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% nexttile
% hold on
% fplot(ddfunMat(1,2,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(ddfunMat(1,2,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(ddfunMat(1,2,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(ddfunMat(1,2,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(ddfunMat(1,2,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% ylabel(char(8706)+""+char(178)+char(8450)+"12/"+char(8706)+char(632)+char(178)+" [GPa]");
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% nexttile
% hold on
% fplot(ddfunMat(3,3,5),[0 1],'Color','#000000','LineWidth',1.5);
% fplot(ddfunMat(3,3,1),[0 1],'-','Color','#D95319','LineWidth',1.5);
% fplot(ddfunMat(3,3,2),[0 1],'--','Color','#D95319','LineWidth',1.5);
% fplot(ddfunMat(3,3,3),[0 1],'-','Color','#0072BD','LineWidth',1.5);
% fplot(ddfunMat(3,3,4),[0 1],'--','Color','#0072BD','LineWidth',1.5);
% ylabel(char(8706)+""+char(178)+char(8450)+"33/"+char(8706)+char(632)+char(178)+" [GPa]");
% xlabel("Damage "+char(632)+" [-]");
% fontsize(gcf,25,'points')
% 
% lg =legend('Analytical','Circle (Volume)','Circle (Length)','Square (Volume)','Square (Length)');
% lg.Layout.Tile = 'East';
% 
%% Force displacement plots (after loading all outputs in output{})
t = tiledlayout(1,2);
nexttile
hold on
plot(output{1}.displacement.value,output{1}.force,'Color','#000000','LineWidth',1.5)
plot(output{2}.displacement.value,output{2}.force,'Color','#D95319','LineWidth',1.5)
plot(output{3}.displacement.value,output{3}.force,'Color','#0072BD','LineWidth',1.5)
xlabel('Displacement [mm]')
xlim([0 0.025])
ylabel('Force [kN]')
ylim([0 2.75])
title('AT1 - Length')

nexttile
hold on
plot(output{4}.displacement.value,output{4}.force,'Color','#000000','LineWidth',1.5)
plot(output{5}.displacement.value,output{5}.force,'Color','#D95319','LineWidth',1.5)
plot(output{6}.displacement.value,output{6}.force,'Color','#0072BD','LineWidth',1.5)
xlabel('Displacement [mm]')
xlim([0 0.04])
ylabel('Force [kN]')
ylim([0 2.75])
title('AT2 - Volume')

fontsize(gcf,25,'points')
lg =legend('Analytical','Circle','Square');
lg.Layout.Tile = 'East';

%% Damage displacement plots
t = tiledlayout(1,2);
nexttile
hold on
plot(output{1}.displacement.value,output{1}.damage.maxValue,'Color','#000000','LineWidth',1.5)
plot(output{2}.displacement.value,output{2}.damage.maxValue,'Color','#D95319','LineWidth',1.5)
plot(output{3}.displacement.value,output{3}.damage.maxValue,'Color','#0072BD','LineWidth',1.5)
xlabel('Displacement [mm]')
xlim([0 0.025])
ylabel("Damage "+char(632)+" [-]")
ylim([0 1])
title('AT1 - Length')

nexttile
hold on
plot(output{4}.displacement.value,output{4}.damage.maxValue,'Color','#000000','LineWidth',1.5)
plot(output{5}.displacement.value,output{5}.damage.maxValue,'Color','#D95319','LineWidth',1.5)
plot(output{6}.displacement.value,output{6}.damage.maxValue,'Color','#0072BD','LineWidth',1.5)
xlabel('Displacement [mm]')
xlim([0 0.04])
ylabel("Damage "+char(632)+" [-]")
ylim([0 1])
title('AT2 - Volume')

fontsize(gcf,25,'points')
lg =legend('Analytical','Circle','Square');
lg.Layout.Tile = 'East';

% %% Damage SEN plots
% t = tiledlayout(2,3);
% nexttile
% output{1}.displacement.value,output{1}.damage.field
% title('Analytical')
% nexttile
% output{2}.displacement.value,output{2}.damage.field
% title('Circle')
% nexttile
% output{3}.displacement.value,output{3}.damage.field
% title('Square')
% nexttile
% output{4}.displacement.value,output{4}.damage.field
% nexttile
% output{5}.displacement.value,output{5}.damage.field
% nexttile
% output{6}.displacement.value,output{6}.damage.field