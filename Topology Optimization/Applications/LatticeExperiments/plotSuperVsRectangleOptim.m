function plotSuperVsRectangleOptim
mCase = 'Middle';
folderR = ['/home/alex/Desktop/ExperimentingPlotRectangle',mCase,'/'];  
[iter,cost,maxStress] = obtainCostShapes(folderR);
f = figure();
h{1} = plot(iter,[cost;maxStress]','b');
hline = findobj(gcf, 'type', 'line');
set(hline(1),'LineStyle','--')

folderS = ['/home/alex/Desktop/ExperimentingPlotSuperEllipse',mCase,'/'];
[iter,costS,maxStressS] = obtainCostShapes(folderS);
hold on
h{2} = plot(iter,[costS;maxStressS]','r');
hline = findobj(gcf, 'type', 'line');
set(hline(1),'LineStyle','--')

leg = {'Rectangle Amplified Stress p', 'Amplified Stress max',...
       'Superellipse Amplified Stress p', 'Superellipse Stress max'};
legend(leg)
printer = plotPrinter(f,h);
path = '/home/alex/Dropbox/GregoireMeetings/GregoireMeeting25Janvier';
output = fullfile(path,['Convergence',mCase]);
printer.print(output)

end





function [iter,stress,maxStress] = obtainCostShapes(folder)
fNameMon = fullfile(folder,'Monitoring.fig');
h = openfig(fNameMon);
handles = findobj(h,'Type','line');
iter = get(handles(5),'Xdata');
stress = get(handles(11),'Ydata');
maxStress = get(handles(4),'Ydata');
close(h)
end