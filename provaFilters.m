file = 'punzon4';
a.fileName = file;
s = FemDataContainer(a);
mesh=s.mesh;
load ("densityIntent5.mat");

f=LagrangianFunction.create(mesh,1,'P1');
f.setFValues(density);

z.mesh=mesh;
z.test=f;
z.trial=f;
% z.filterType = 'PDE';
% filter = Filter.create(z);
% fFilter= filter.compute(f,2);
% fFilter.print('hola');


zMax = max(mesh.coord(:,3));

% bolt domain
c1x= 584.37;    c1y= 28.562;
c2x= 584.37;    c2y= 68.562;
hBolt=60;       rBolt=6.5;

bolt1 = @(x) (x(:,1)-c1x).^2 + (x(:,2)-c1y).^2 <= rBolt^2 & ...
    x(:,3) >= hBolt & x(:,3) <= zMax;
bolt2 = @(x) (x(:,1)-c2x).^2 + (x(:,2)-c2y).^2 <= rBolt^2 & ...
    x(:,3) >= hBolt & x(:,3) <= zMax;

% Guides and bottom fixed
isBottom = @(x) x(:,3)<= -16;
% guide1   = @(x) x(:,2)<= yMin+10.5;
% guide2   = @(x) x(:,2)>= yMax-10.5;
guide1   = @(x) x(:,2)<= 20.179;
guide2   = @(x) x(:,2)>= 76.729;
isFixed  = computeFixedVolumeDomain(mesh,...
    @(x) guide1(x) | guide2(x) | isBottom(x)| bolt1(x) | bolt2(x),...
    'Density');


isFixVals = zeros(size(f.fValues));
isFixVals(isFixed.nodes) = isFixed.values;
z.mesh=mesh;
z.test=f;
z.trial=f;
z.filterType = 'PDEDir';
z.isFixed = LagrangianFunction.create(mesh,1,'P1');
z.isFixed.setFValues(isFixVals);
filter = Filter.create(z);
fFilter2= filter.compute(f,2);
fFilter2.print('PDEsol');






function isFixed = computeFixedVolumeDomain(mesh,cond,type)
    coor  = mesh.coord;
    nodes = find(cond(coor));
    isFixed.nodes = nodes;
    switch type
        case 'Density'
            values = ones(size(nodes));
        case 'LevelSet'
            values = -ones(size(nodes));
    end
    isFixed.values = values;
end