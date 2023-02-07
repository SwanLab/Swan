function ComputingTopOpt_ComplianceStudy_AL
clear
clc
close all

% set(0,'DefaultFigureVisible','off');

penaltyVec = [1];
trustVec = [1]; % 0.1

s.testName = 'PerimeterAsConstraint';
for i = 1:length(penaltyVec)
    for j = 1:length(trustVec)
        close all
        k.penaltyProv = penaltyVec(i);
        k.trustProv = trustVec(j);
        t = TopOptComputer(s);
        t.compute(k);
        trustName = num2str(trustVec(j));
        trustName(trustName=='.')='d';
        namerho = ['rho',num2str(penaltyVec(i))];
        nametrust{i,j} = ['trust',trustName];
        name{i,j} = [namerho,nametrust{i,j}];
        saveas(gcf,['ResultsInkscape/',name{i,j},'.fig'])
    end
end

% s.testName = 'test_cantilever2';%''testJose';
% s.testName = 'test_cantilever_nullspace';
% s.testName = 'PerimeterAsConstraint';
% k.penaltyProv = 40;
% k.trustProv   = 100;
% t = TopOptComputer(s);
% t.compute(k);

for i = 1:length(penaltyVec)
    for j = 1:length(trustVec)
        close all
        fig = openfig(['ResultsInkscape/',name{i,j},'.fig']);
        dataObjs = findobj(fig,'-property','YData');
        nPlots = length(fig.Children);
        for p = 1:nPlots
            switch fig.Children(p).Title.String
                case 'Compliance non scaled'
                    YVolum{1,:,i,j} = dataObjs(p).YData;
                case 'Line Search'
                    Ytau{1,:,i,j} = dataObjs(p).YData;
                case '\lambda_V_o_l_u_m_e_C_o_n_s_t_r_a_i_n_t'
                    Ylambda{1,:,i,j} = dataObjs(p).YData;
                case 'Cost'
                    YCompliance{1,:,i,j} = dataObjs(p).YData;
            end
        end
        XIter{1,:,i,j} = 1:1:length(YVolum{1,:,i,j});
    end
end
close all

for j = 1:length(trustVec)
    figure
    for i = 1:length(penaltyVec)
        plot(XIter{1,:,i,j},YVolum{1,:,i,j})
        hold on
    end
    xlabel('Number of iterations','Interpreter','latex')
    ylabel('Volume','Interpreter','latex')
    legend('$\rho=1$','$\rho=10$','$\rho=100$','Interpreter','latex')
    grid on
    grid minor
    saveName = ['ResultsInkscape/Postprocess/Volum',nametrust{i,j},'.svg'];
    saveas(gcf,saveName)
    close all
end

for j = 1:length(trustVec)
    figure
    for i = 1:length(penaltyVec)
        plot(XIter{1,:,i,j},Ytau{1,:,i,j})
        hold on
    end
    xlabel('Number of iterations','Interpreter','latex')
    ylabel('Line search','Interpreter','latex')
    legend('$\rho=1$','$\rho=10$','$\rho=100$','Interpreter','latex')
    grid on
    grid minor
    saveName = ['ResultsInkscape/Postprocess/LineSearch',nametrust{i,j},'.svg'];
    saveas(gcf,saveName)
    close all
end

for j = 1:length(trustVec)
    figure
    for i = 1:length(penaltyVec)
        plot(XIter{1,:,i,j},Ylambda{1,:,i,j})
        hold on
    end
    xlabel('Number of iterations','Interpreter','latex')
    ylabel('$\lambda$','Interpreter','latex')
    legend('$\rho=1$','$\rho=10$','$\rho=100$','Interpreter','latex')
    grid on
    grid minor
    saveName = ['ResultsInkscape/Postprocess/Lambda',nametrust{i,j},'.svg'];
    saveas(gcf,saveName)
    close all
end

for j = 1:length(trustVec)
    figure
    for i = 1:length(penaltyVec)
        plot(XIter{1,:,i,j},YCompliance{1,:,i,j})
        hold on
    end
    xlabel('Number of iterations','Interpreter','latex')
    ylabel('Compliance','Interpreter','latex')
    legend('$\rho=1$','$\rho=10$','$\rho=100$','Interpreter','latex')
    grid on
    grid minor
    saveName = ['ResultsInkscape/Postprocess/Perimeter',nametrust{i,j},'.svg'];
    saveas(gcf,saveName)
    close all
end

set(0,'DefaultFigureVisible','on');

end