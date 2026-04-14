function obtainGIF(gifName,designVariable,nIter)
set(0,'DefaultFigureVisible','off');
deltaTime = 0.01;
m = designVariable.fun.mesh;
xmin = min(m.coord(:,1));
xmax = max(m.coord(:,1));
ymin = min(m.coord(:,2));
ymax = max(m.coord(:,2));

f = designVariable.fun.fValues;
switch designVariable.type
    case 'LevelSet'
        uMesh = designVariable.getUnfittedMesh();
        uMesh.compute(f);
        gifFig = figure;
        uMesh.plotStructureInColor('black');
    case 'Density'
        p1.mesh    = m;
        p1.fValues = f;
        p1.order   = 'P1';
        RhoNodal   = LagrangianFunction(p1);
        q = Quadrature.create(m,0);
        xV = q.posgp;
        RhoElem = squeeze(RhoNodal.evaluate(xV));

        gifFig = figure;
        axis off
        axis equal
        axes = gifFig.Children;
        patchHandle = patch(axes,'Faces',m.connec,'Vertices',m.coord,...
            'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
        set(axes,'ALim',[0, 1],'XTick',[],'YTick',[]);
        set(patchHandle,'FaceVertexAlphaData',RhoElem,'FaceAlpha','flat');
end
hold on
fig = gifFig;
fig.CurrentAxes.XLim = [xmin xmax];
fig.CurrentAxes.YLim = [ymin ymax];
axis([xmin xmax ymin ymax])
gifname = [gifName,'.gif'];
set(gca, 'Visible', 'off')

frame = getframe(fig);
[A,map] = rgb2ind(frame.cdata,256);
if nIter == 0
    imwrite(A,map,gifname,"gif","LoopCount",0,"DelayTime",deltaTime);
else
    imwrite(A,map,gifname,"gif","WriteMode","append","DelayTime",deltaTime);
end
close(gifFig);
set(0,'DefaultFigureVisible','on');
end