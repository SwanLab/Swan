clear

% Cube chair BCs
run('TFGAlvaro/GiD/CubeChairAlvaro.gid/CubeChairAlvaro.m')

lnodes             = [];
pointload_complete = [];
coor               = gidcoord(:,2:end);

% Nodes constrained
t                 = 0.2; % Hole size
x                 = coor(:,1);
y                 = coor(:,2);
z                 = coor(:,3);
xMax              = max(x);
yMax              = max(y);
nConstz           = find(z==0);
nConsty           = find(y<=t*yMax | y>=(1-t)*yMax);
nConstx           = find(x<=t*xMax | x>=(1-t)*xMax);
nConstyz          = nConstz(ismember(nConstz,nConsty,'rows'));
nConstxz          = nConstz(ismember(nConstz,nConstx,'rows'));
nConst            = nConstyz(ismember(nConstyz,nConstxz,'rows'));
nConstRepeated    = [nConst';nConst';nConst'];
lnodes            = [lnodes,nConstRepeated(:)];
lnodes(:,[2,3])   = 0;
lnodes(1:3:end,2) = 1;
lnodes(2:3:end,2) = 2;
lnodes(3:3:end,2) = 3;

% External forces
t         = 0.3; % Hole size (0.3A x 0.3A)
zTip      = max(z);
ref       = min(yMax,xMax);
nForcesz  = find(z==zTip);
nForcesy  = find(y>=(yMax-t*ref)/2 & y<=(yMax+t*ref)/2);
nForcesx  = find(x>=(xMax-t*ref)/2 & x<=(xMax+t*ref)/2);
nForcesyz = nForcesz(ismember(nForcesz,nForcesy,'rows'));
nForcesxz = nForcesz(ismember(nForcesz,nForcesx,'rows'));
nForces   = nForcesyz(ismember(nForcesyz,nForcesxz,'rows'));

nForcesRepeated               = [nForces';nForces';nForces'];
pointload_complete            = [pointload_complete,nForcesRepeated(:)];
pointload_complete(:,[2,3])   = 0;
pointload_complete(1:3:end,2) = 1;
pointload_complete(2:3:end,2) = 2;
pointload_complete(3:3:end,2) = 3;
pointload_complete(3:3:end,3) = -1;

save('TFGAlvaro/BoundaryConditions/dataCubeChair.mat','lnodes','pointload_complete');