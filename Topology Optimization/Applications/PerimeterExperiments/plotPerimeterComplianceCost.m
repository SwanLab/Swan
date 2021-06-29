function plotPerimeterComplianceCost
%perCase = 'Total';
perCase = 'Relative';

plotCost(perCase);
end

function plotCost(perCase)
folder = ['/media/alex/My Passport/PerimeterResults/FineMesh/',perCase,'Perimeter01/'];  

[iter,compliance,perimeter] = obtainCostShapes(folder);

f = figure();
hold on
pN{1} = plot(iter,0.1*perimeter);
pN{2} = plot(iter,0.4*ones(size(iter)));
yyaxis right
pN{3} = plot(iter,compliance);

leg = legend({['$\textrm{',perCase,' perimeter } \alpha \textrm{Per}^{T}_{\varepsilon}(\Omega)$'],'$\textrm{Volume }\textrm{V}(\Omega)$','$ \textrm{Compliance } \textrm{C}(\Omega) $'});
set(leg,'Interpreter','latex')
p = plotPrinter(f,pN);
p.print(fullfile(folder,['Cost']));

end

function [iter,compliance,perimeter] = obtainCostShapes(folder)
fNameMon = fullfile(folder,'Monitoring.fig');
h = openfig(fNameMon);
handles = findobj(h,'Type','line');
iter = get(handles(5),'Xdata');
compliance = get(handles(6),'Ydata');
perimeter = get(handles(5),'Ydata');
close all
end