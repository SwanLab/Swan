%% Goal: delete integrators
% Previous goal: delete BMatrixComputer
clc; clear; close all;

% Test and trial
mesh = UnitTriangleMesh(10,10);
quad = Quadrature.set(mesh.type);
quad.computeQuadrature('QUADRATIC');
xV = quad.posgp;

ndimf = 2;
trial = TrialFunction.create(mesh, ndimf, 'P1');
test  = TestFunction.create(mesh, ndimf, 'P1');

%% Mass Matrix
shTr = Shape(trial).evaluate(xV);
shTe = Shape(test).evaluate(xV);
dVolu(1,:,:)  = mesh.computeDvolume(quad);
mult = pagemtimes(shTr',shTe);
integ = bsxfun(@times, mult, dVolu);
% kron(squeeze(A(:,:,1)),eye(2)) % !!!!
integ = openprod(integ,eye(ndimf));

M = massmatrix(quad,mesh,test,trial);

shDerTe = ShapeDer(test).evaluate(xV);

%% Fake stiffness matrix

shTr = ShapeDerSym(trial).evaluate(xV);
shTe = ShapeDerSym(test).evaluate(xV);
dVolu(1,:,:)  = mesh.computeDvolume(quad);
mult = pagemtimes(permute(shTr,[2 1 3 4]),shTe);


%% Functions
function z = openprod(A,B)
    z = bsxfun(@times, permute(A, [4 1 5 2 3]), permute(B, [1 4 2 5 3]));  %// step 1
    z = reshape(z, size(A,1)*size(B,1), size(A,2)*size(B,2), size(A,3)); 
end

function lhs = massmatrix(quad,mesh, test, trial)
    xV   = quad.posgp;
    shapesTest  = test.computeShapeFunctions(xV);
    shapesTrial = trial.computeShapeFunctions(xV);
    dVolu  = mesh.computeDvolume(quad);

    nGaus  = quad.ngaus;
    nElem  = size(dVolu,2);
    nNodeTest  = size(shapesTest,1);
    nNodeTrial = size(shapesTrial,1);
    nDofTest   = nNodeTest*test.ndimf;
    nDofTrial  = nNodeTrial*trial.ndimf;

    M = zeros(nDofTest, nDofTrial, nElem);
    for igauss = 1 :nGaus
        for inode= 1:nNodeTest
            for jnode= 1:nNodeTrial
                for iunkn= 1:test.ndimf
               %     for junkn= 1:obj.trial.ndimf
                        idof = test.ndimf*(inode-1)+iunkn;
                        jdof = trial.ndimf*(jnode-1)+iunkn;
                        dvol = dVolu(igauss,:);
                        Ni = shapesTest(inode,igauss,:);
                        Nj = shapesTrial(jnode,igauss,:);
                        v = squeeze(Ni.*Nj);
                        M(idof, jdof, :)= squeeze(M(idof,jdof,:)) ...
                            + v(:).*dvol';
               %     end
                end
            end
        end
    end
    lhs = M;
end


function lhs = stiffnessmatrix(quad,mesh, test, trial)
    xV = quad.posgp;
    dNdxTs = test.evaluateCartesianDerivatives(xV);
    dNdxTr = trial.evaluateCartesianDerivatives(xV);
    dVolu = mesh.computeDvolume(obj.quadrature);
    nGaus = quad.ngaus;
    nElem = size(dVolu,2);

    nNodETs = size(dNdxTs,2);
    nDofETs = nNodETs*test.ndimf;
    nNodETr = size(dNdxTr,2);
    nDofETr = nNodETr*trial.ndimf;

    BcompTs = createBComputer(test, dNdxTs);
    BcompTr = createBComputer(trial, dNdxTr);
    lhs = zeros(nDofETs,nDofETr,nElem);
    for igaus = 1:nGaus
        BmatTs = BcompTs.compute(igaus);
        BmatTr = BcompTr.compute(igaus);
        dV(1,1,:) = dVolu(igaus,:)';
        Bt   = permute(BmatTs,[2 1 3]);
        BtCB = pagemtimes(Bt, BmatTr);
        lhs = lhs + bsxfun(@times, BtCB, dV);
    end
end

function Bcomp = createBComputer( fun, dNdx)
    s.fun  = fun;
    s.dNdx = dNdx;
    Bcomp = BMatrixComputer(s);
end