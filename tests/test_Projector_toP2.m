sM.coord = [0 0;1 0;0 1];
sM.connec = [1 2 3];

m = Mesh(sM);
mesh = m;

sAF.fHandle = @(x) x(1,:,:); % f(x) = x
sAF.ndimf   = 1; % number of dimensions
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);




p2fun = xFun.project('P2');
p2fun.plot()


xFun.evualaute(xG);
p2.evaluate(xG);

for iH = 1:nH
for iOrder = 1:2
    
        createMesh()


(xFun - p2).^2

end


%ToDo
%1. Merge P1 and P2 by LagrangianFeFunction and check tests (AllTests)
%using Dofs name
%2. Convergence analysis for P1 and P2 in Tutorial