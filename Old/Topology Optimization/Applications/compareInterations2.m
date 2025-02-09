function compareInterations2

fCase{1}  = 'ExperimentingPlot';
folder{1} = '/media/alex/My Passport/LatticeResults/StressNormSuperEllipseRotation';

fCase{2} = 'ExperimentingPlot';
folder{2} = '/media/alex/My Passport/LatticeResults/StressNormRectangleRotation';


%fCase{3} = 'LatticeExperimentInputCantileverSymmetricMeshSuperEllipsePDE';
%folder{3} = '/media/alex/My Passport/LatticeResults/CantileverSymmetricMeshSuperEllipsePDE10000iteration';

f = figure();
hold on
it = 1;
for iPlot = 1:numel(fCase)
    fName = fullfile(folder{iPlot},fCase{iPlot});
    [xV,yV,cV] = getCostFunctions(folder{iPlot});
   % for iline = 1:size(yV,2)
        p{it} = plot(xV,yV(:,5));
        it = it + 1;
  %  end
end
legend('SuperEllipse','Rectangle')
pr = plotPrinter(f,p);
pr.print(fullfile('/home/alex/Dropbox/GregMeeting30Octubre',['Iterations3']));
end


function [xV,yV,cV] = getCostFunctions(fName)
fNameMon = fullfile([fName,'/Monitoring.fig']);
h = openfig(fNameMon);
handles = findobj(h,'Type','line');
xV = get(handles(1),'Xdata');
for i = 1:numel(handles)
    yV(:,i) = get(handles(i),'Ydata');
    cV(:,i) = get(handles(i),'Color') ;
end
close(h)
end