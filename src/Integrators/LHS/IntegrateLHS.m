function LHS = IntegrateLHS(f,test,trial,mesh,quadOrder)
if nargin < 5 || isempty(quadOrder)
    qTe = test.getOrderNum();
    qTr = trial.getOrderNum();
    quadOrder = qTe + qTr;
end
lhs = integrateElementalLHS(f,test,trial,mesh,quadOrder);
LHS = assembleMatrix(lhs,test, trial);
end

function lhs = integrateElementalLHS(f,test,trial,mesh,quadOrder)
quad = Quadrature.create(mesh,quadOrder);
xV = quad.posgp;
w  = quad.weigp;
lhs    = zeros(test.nDofsElem,trial.nDofsElem,mesh.nelem);
detJ   = DetJacobian(mesh);
v = @(i) Test(test,i);
u = @(j) Test(trial,j);
for i = 1:test.nDofsElem
    for j = 1:trial.nDofsElem
        int = (f(u(j),v(i)).*detJ)*w';
        lhs(i,j,:) = lhs(i,j,:) + int.evaluate(xV);
    end
end
end

