addpath(genpath(pwd))

s.geometryType = "Surface";
% s.coord = [0,0;1,0;0,1;1,1];
% s.connec = [1 2 3;
%             2 3 4];
s.coord = [0,0;1,0;0,1];
s.connec = [1 2 3];
m = Mesh.create(s);

rt = RaviartThomasFunction.create(m,1,1);

sAF.fHandle = @(x) [x(1,:,:) ; x(2,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = m;
xFun = AnalyticalFunction(sAF);

% p1fun = xFun.project('P1');
rtfun = xFun.project('RT');