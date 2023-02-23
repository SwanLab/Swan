function tryingLHS
n = 50;
xCoord =  sort(rand(1,n));
s.coord(:,1) = xCoord;
s.coord(:,2) = zeros(n,1);
s.connec(:,1) = 1:n-1;
s.connec(:,2) = 2:n;
s.kFace = -1;
m = Mesh(s);
m.plot();

s.mesh               = m;
s.ndimf              = 1;
s.interpolationOrder = 'LINEAR';
f = Field(s);

l.mesh = m;
l.field = f;
l.type = 'StiffnessMatrix';
lhs = LHSintegrator.create(l);
K = lhs.compute();

s.type  = 'MassMatrix';
s.mesh  = m;
s.field = f;
LHS     = LHSintegrator.create(s);
M   = LHS.compute();

e = 0.1*m.computeMeanCellSize();
LHS = e*K + M;

xmean = (max(xCoord)-min(xCoord))/2;
f = heaviside(xCoord-xmean);

s.mesh = m;
s.fValues = f;
fL2 = P1Function(s);
fL2.plot
RHS = M*(fL2.fValues)';


s.type    =  'DIRECT';
solver    = Solver.create(s);
fH1Values = solver.solve(LHS,RHS);

s.mesh = m;
s.fValues = fH1Values;
fH1 = P1Function(s);
fH1.plot()

m.plot();
hold on
plot(xCoord,fL2.fValues,'-+');
plot(xCoord,fH1.fValues,'-+');
end