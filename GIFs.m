clear
set(0,'DefaultFigureVisible','off');

% INPUTS
s.filename = 'romeCase222';
testname = 'bridgeRome';
s.designVariable = 'Density';
s.initial = 528897;
s.final = 558547;
stepVec = [33:895,1000:1095,1200:1280,1400:1485,1600:1802];
% END INPUTS

m.fileName = testname;
m = FemDataContainer(m);
m = m.mesh;
s.mesh  = m;

switch s.designVariable
    case 'LevelSet'
        sB.backgroundMesh = m;
        sB.dimension = 1:3;
        sB.type = 'FromReactangularBox';
        bMc = BoundaryMeshCreator.create(sB);
        boundaryMesh  = bMc.create();
        sU.backgroundMesh = m;
        sU.boundaryMesh   = boundaryMesh;
        uMesh = UnfittedMesh(sU);
        s.uMesh = uMesh;
    case 'Density'
        % UnfittedMesh not required
end

count = 1;
for j=1:length(stepVec)
    s.count = count;
    s.step = stepVec(j);
    snapshot(s);
    count = count+1;
end

set(0,'DefaultFigureVisible','on');

function snapshot(cParams)

step     = char(string(cParams.step));
filename = [cParams.filename,step,'.flavia.res'];
initial  = cParams.initial;
final    = cParams.final;

file = fopen(filename);
a    = textscan(file,'%s','delimiter','\n');
a    = a{1};

f = [];

for i=initial:final
    b = a{i,1};
    c = not(b==' ');
    d = find(c==false);
    e = d(1)+1;
    f = [f;str2double(b(e:end))];
end

m = cParams.mesh;
xmin = min(m.coord(:,1));
xmax = max(m.coord(:,1));
ymin = min(m.coord(:,2));
ymax = max(m.coord(:,2));

switch cParams.designVariable
    case 'LevelSet'
        uMesh = cParams.uMesh;
        uMesh.compute(f);
        uMesh.plotStructureInColor('black');
        hold on
    case 'Density'
        p1.connec  = m.connec;
        p1.type    = m.type;
        p1.fValues = f;
        RhoNodal   = P1Function(p1);
        p1.mesh    = m;
        toP0       = Projector_toP0(p1);
        RhoElemFun = toP0.project(RhoNodal);
        RhoElem    = squeeze(RhoElemFun.fValues);

        figHandle  = figure;
        axis off
        axis equal
        axes = figHandle.Children;
        patchHandle = patch(axes,'Faces',m.connec,'Vertices',m.coord,...
            'EdgeColor','none','LineStyle','none','FaceLighting','none' ,'AmbientStrength', .75);
        set(axes,'ALim',[0, 1],'XTick',[],'YTick',[]);
        set(patchHandle,'FaceVertexAlphaData',RhoElem,'FaceAlpha','flat');
end

fig = gcf;
fig.CurrentAxes.XLim = [xmin xmax];
fig.CurrentAxes.YLim = [ymin ymax];
axis([xmin xmax ymin ymax])
gifname = [cParams.filename,'.gif'];
set(gca, 'Visible', 'off')

frame = getframe(fig);
[A,map] = rgb2ind(frame.cdata,256);
if cParams.count == 1
    imwrite(A,map,gifname,"gif","LoopCount",0,"DelayTime",0.006);
else
    imwrite(A,map,gifname,"gif","WriteMode","append","DelayTime",0.006);
end
close all

end