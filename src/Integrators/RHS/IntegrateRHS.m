function RHS = IntegrateRHS(f,test,mesh,quadOrder)
if nargin < 4 || isempty(quadOrder)
    quadOrder = 2;
end
quad = Quadrature.create(mesh,quadOrder);
xV = quad.posgp;
w  = quad.weigp;
rhs    = zeros(test.nDofsElem,mesh.nelem);
detJ   = DetJacobian(omesh);
v = @(i) Test(test,i);
for i = 1:test.nDofsElem
    int = (f(v(i)).*detJ)*w';
    rhs(i,:) = rhs(i,:) + squeezeParticular(int.evaluate(xV),2);
end
RHS = assembleVector(rhs, test);
end