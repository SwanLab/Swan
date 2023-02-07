function ComputingTopOpt_ComplianceStudy_NS
clear
clc
close all

% set(0,'DefaultFigureVisible','off');

fVec = 1;

s.testName = 'PerimeterAsConstraint';
for i = 1:length(fVec)
    close all
    k.fProv = fVec(i);
    t = TopOptComputer(s);
    t.compute(k);
    fName = num2str(fVec(i));
    fName(fName=='.')='d';
    namef = ['f',fName];
    name{i} = [namef];
    saveas(gcf,['ResultsInkscape/',name{i},'.fig'])
end

for i = 1:length(fVec)
    close all
    fig = openfig(['ResultsInkscape/',name{i},'.fig']);
    dataObjs = findobj(fig,'-property','YData');
    nPlots = length(fig.Children);
    for p = 1:nPlots
        switch fig.Children(p).Title.String
            case 'Compliance non scaled'
                YVolum{1,:,i} = dataObjs(p).YData;
            case 'Line Search'
                Ytau{1,:,i} = dataObjs(p).YData;
            case '\lambda_V_o_l_u_m_e_C_o_n_s_t_r_a_i_n_t'
                Ylambda{1,:,i} = dataObjs(p).YData;
            case 'Cost'
                YCompliance{1,:,i} = dataObjs(p).YData;
        end
    end
    XIter{1,:,i} = 1:1:length(YVolum{1,:,i});
end
close all

figure
for i = 1:length(fVec)
    plot(XIter{1,:,i},YVolum{1,:,i})
    hold on
end
xlabel('Number of iterations','Interpreter','latex')
ylabel('Volume','Interpreter','latex')
legend('$f=0.01$','$f=0.1$','$f=1$','$f=10$','$f=100$','Interpreter','latex')
grid on
grid minor
saveName = ['ResultsInkscape/Postprocess/Volum','.svg'];
saveas(gcf,saveName)
close all

figure
for i = 1:length(fVec)
    plot(XIter{1,:,i},Ytau{1,:,i})
    hold on
end
xlabel('Number of iterations','Interpreter','latex')
ylabel('Line search','Interpreter','latex')
legend('$f=0.01$','$f=0.1$','$f=1$','$f=10$','$f=100$','Interpreter','latex')
grid on
grid minor
saveName = ['ResultsInkscape/Postprocess/LineSearch','.svg'];
saveas(gcf,saveName)
close all

figure
for i = 1:length(fVec)
    plot(XIter{1,:,i},Ylambda{1,:,i})
    hold on
end
xlabel('Number of iterations','Interpreter','latex')
ylabel('$\lambda$','Interpreter','latex')
legend('$f=0.01$','$f=0.1$','$f=1$','$f=10$','$f=100$','Interpreter','latex')
grid on
grid minor
saveName = ['ResultsInkscape/Postprocess/Lambda','.svg'];
saveas(gcf,saveName)
close all

figure
for i = 1:length(fVec)
    plot(XIter{1,:,i},YCompliance{1,:,i})
    hold on
end
xlabel('Number of iterations','Interpreter','latex')
ylabel('Compliance','Interpreter','latex')
legend('$f=0.01$','$f=0.1$','$f=1$','$f=10$','$f=100$','Interpreter','latex')
grid on
grid minor
saveName = ['ResultsInkscape/Postprocess/Perimeter','.svg'];
saveas(gcf,saveName)
close all

set(0,'DefaultFigureVisible','on');

end