clear

% Cube chair BCs
run('Benchmarking/ChairAlvaro.m')

lnodes             = [];
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

% External forces on horizontal surface
fHS       = 3;
zHS       = 0.5*max(z);
nForcesz  = find(z==zHS);
nForcesx  = find(x<=xMax-8);
nForces   = nForcesz(ismember(nForcesz,nForcesx,'rows'));

pointload_complete1            = [];
nForcesRepeated                = [nForces';nForces';nForces'];
pointload_complete1            = [pointload_complete1,nForcesRepeated(:)];
pointload_complete1(:,[2,3])   = 0;
pointload_complete1(1:3:end,2) = 1;
pointload_complete1(2:3:end,2) = 2;
pointload_complete1(3:3:end,2) = 3;
pointload_complete1(3:3:end,3) = -fHS;


% External forces on vertical surface
fVS = 1;
xVS = xMax-8;
nForcesx = find(x==xVS);
nForcesz = find(z>=zHS);
nForces  = nForcesx(ismember(nForcesx,nForcesz,'rows'));

pointload_complete2            = [];
nForcesRepeated                = [nForces';nForces';nForces'];
pointload_complete2            = [pointload_complete2,nForcesRepeated(:)];
pointload_complete2(:,[2,3])   = 0;
pointload_complete2(1:3:end,2) = 1;
pointload_complete2(2:3:end,2) = 2;
pointload_complete2(3:3:end,2) = 3;
pointload_complete2(1:3:end,3) = fVS;

pointload_complete = [pointload_complete1;pointload_complete2];


save('Benchmarking/dataChair.mat','lnodes','pointload_complete');