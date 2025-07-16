function ObtainGif(gifName,variable)
deltaTime = 0.01;
m         = variable.fun.mesh;
xmin      = min(m.coord(:,1));
xmax      = max(m.coord(:,1));
ymin      = min(m.coord(:,2));
ymax      = max(m.coord(:,2));
fun       = variable.fun;
q = Quadrature.create(m,0);
xV = q.posgp;
RhoElem = squeeze(fun.evaluate(xV));
gifFig = figure();
axis off
axis equal
axes = gifFig.Children;
patchHandle = patch(axes,'Faces',m.connec,'Vertices',m.coord,...
    'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
set(axes,'ALim',[0, 1],'XTick',[],'YTick',[]);
set(patchHandle,'FaceVertexAlphaData',RhoElem,'FaceAlpha','flat');

hold on;
fig = gifFig;
fig.CurrentAxes.XLim = [xmin xmax];
fig.CurrentAxes.YLim = [ymin ymax];
axis([xmin xmax ymin ymax])
gifname = [gifName,'.gif'];
set(gca, 'Visible', 'off')

frame = getframe(fig);
[A,map] = rgb2ind(frame.cdata,256);
if obj.nIter == 0
    imwrite(A,map,gifname,"gif","LoopCount",0,"DelayTime",deltaTime);
else
    imwrite(A,map,gifname,"gif","WriteMode","append","DelayTime",deltaTime);
end
close(gifFig);
end