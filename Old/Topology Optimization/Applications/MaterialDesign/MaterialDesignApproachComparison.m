function MaterialDesignApproachComparison

matCase = 'Horizontal';
path = '/media/alex/My Passport/MaterialDesign/';
element = {'Tri','Quad'};
desVar  = {'Density','LevelSet'};
filter  = {'P1';'PDE'};

fCase{1} = 'CompositeMaterialDesign';
folder{1} = fullfile(path,matCase,element,desVar,filter);


fCase{2} = 'CompositeMaterialDesign';
folder{2} = '/media/alex/My Passport/MaterialDesign/Quadrilater/Horizontal/Density';

%fCase{2} = 'ExperimentingPlot';
%folder{2} = '/media/alex/My Passport/LatticeResults/StressNormRectangleRotation';

f = figure();
hold on
it = 1;
for iPlot = 1:numel(fCase)
   % fName = fullfile(folder{iPlot},fCase{iPlot});
    [xV,yV,cV] = getCostFunctions(folder{iPlot});
   % for iline = 1:size(yV,2)
        p{it} = plot(xV,yV(:,7));
        it = it + 1;
  %  end
end
legend('LevelSet','Density')
pr = plotPrinter(f,p);
pr.print(fullfile('/home/alex/Dropbox/MaterialDesign',['HorizontalQuadrilater']));
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