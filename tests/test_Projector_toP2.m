sM.coord = [0 0;1 0;0 1];
sM.connec = [1 2 3];

m = Mesh(sM);
mesh = m;

sAF.fHandle = @(x) x(1,:,:); % f(x) = x
sAF.ndimf   = 1; % number of dimensions
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p2fun = xFun.project('P2');