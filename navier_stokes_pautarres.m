mesh = UnitQuadMesh(10,10);

u = LagrangianFunction.create(mesh,2,'P2');
p = LagrangianFunction.create(mesh,2,'P1');

isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

sDir{1}.domain    = @(coor) isLeft(coor) | isRight(coor) | isTop(coor) | isBottom(coor);
sDir{1}.direction = [1,2];
sDir{1}.value     = 0;

sDir{2}.domain    = @(coor) isRight(coor) & isTop(coor);
sDir{2}.direction = 1;
sDir{2}.value     = 0;

velocityBC.domain = 'Border';
velocityBC.value  = 0;

inBC.pressure  = cParams.bc.pressure;
inBC.velocity  = cParams.bc.velocity;
inBC.pointload = [];
inBC.velocityBC    = cParams.bc.velocityBC;
inBC.forcesFormula = cParams.bc.forcesFormula;
obj.inputBC    = inBC;






s1.type  = 'StiffnessMatrix';
s1.mesh  = mesh;
s1.test  = u;
s1.trial = u;
LHS1 = LHSintegrator.create(s1);
K = LHS1.compute();

s2.type = 'WeakDivergence';
s2.mesh = mesh;
s2.trial = p;
s2.test  = u;
LHS2 = LHSintegrator.create(s2);
D = LHS2.compute();

sz = size(D, 2);
BB = sparse(sz,sz);

LHS = [K D; D', BB];



a.type = 'ShapeFunction';
a.mesh = mesh;
a.quadType = 3;
rhsI = RHSintegrator.create(a);
test = u;
Fext = rhsI.compute(obj.forcesFormula,test);

