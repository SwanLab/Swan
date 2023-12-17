clear

% Cantilever BCs
run('TFGAlvaro/GiD/CantileverAlvaro.gid/CantileverAlvaro.m')

lnodes             = [];
pointload_complete = [];
coor               = gidcoord(:,2:end);

% Nodes constrained
x                 = coor(:,1);
nConst            = find(x==0);
nConstRepeated    = [nConst';nConst';nConst'];
lnodes            = [lnodes,nConstRepeated(:)];
lnodes(:,[2,3])   = 0;
lnodes(1:3:end,2) = 1;
lnodes(2:3:end,2) = 2;
lnodes(3:3:end,2) = 3;

% External forces
t2        = 0.6;
t         = 0.3; % Hole size (0.3A x 0.3A)
xTip      = max(x);
y         = coor(:,2);
z         = coor(:,3);
yMax      = max(y);
zMax      = max(z);
ref       = min(yMax,zMax);
nForcesx  = find(x==xTip);
nForcesy  = find(y>=(yMax-t*ref)/2 & y<=(yMax+t*ref)/2);
nForcesz  = find(z>=(zMax-t2*ref)/2 & z<=(zMax+t2*ref)/2);
nForcesyx = nForcesx(ismember(nForcesx,nForcesy,'rows'));
nForceszx = nForcesx(ismember(nForcesx,nForcesz,'rows'));
nForces   = nForcesyx(ismember(nForcesyx,nForceszx,'rows'));

nForcesRepeated               = [nForces';nForces';nForces'];
pointload_complete            = [pointload_complete,nForcesRepeated(:)];
pointload_complete(:,[2,3])   = 0;
pointload_complete(1:3:end,2) = 1;
pointload_complete(2:3:end,2) = 2;
pointload_complete(3:3:end,2) = 3;
pointload_complete(3:3:end,3) = -1;

save('TFGAlvaro/BoundaryConditions/dataCantilever.mat','lnodes','pointload_complete');